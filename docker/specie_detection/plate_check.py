#!/usr/bin/env python3
"""plate_check.py - compare sequencing results with the sending file, plate by plate.

The sending file (xlsx) lists, before sequencing, which sample sits in which plate/well and
which species it should be. After sequencing, the species detected for each sample
(kraken2 or ML) is compared with that expectation, to spot sample mix-ups quickly.

Output: one xlsx with one sheet per plate THAT HAS RESULTS (an 8x12 grid of sample names):
    green  = detected species matches the sending file
    orange = same genus, different species (approximate match)
    red    = does not match (also listed under the grid: sample, detected, expected, well)
    blue   = unsure: a known confusion of the scanner (e.g. Shigella detected as E. coli, see CONFUSED_PAIRS)
    grey   = no result for this sample (not run), or no expected species in the sending file
plus a "Notes" sheet (summary, normalised species names, warnings).

Standard library only: nothing to pip install.

Usage:
    python plate_check.py --sending sending_file.xlsx --results results.tsv --out plate_check.xlsx

--results is a TSV/CSV exported from Terra (or any 2-column file): an id column and a column with
the detected species (auto-detected: "...species_detected"); use --id-col / --species-col to force.
"""
import argparse
import csv
import difflib
import io
import re
import sys
import unicodedata
import zipfile
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape

# ====================================================================================
# Settings you may want to edit
# ====================================================================================

# Accepted header names in the sending file (compared in lower case, without accents or trailing dots).
# The table can start on any row: the first row that contains all four columns is the header.
HEADER_ALIASES = {
    "code": ["samples codes", "sample codes", "sample code", "samples code", "sample id",
             "sample name", "sample", "code"],
    "species": ["species", "espece", "organism"],
    "plate": ["plate no", "plate number", "plate"],
    "well": ["well no", "well number", "well", "position"],
}

# Species used to expand abbreviations ("S. aureus") and to fix typos ("pertusiss").
# Add yours here. Detected names found in the results file are added automatically.
KNOWN_SPECIES = [
    "Salmonella enterica", "Listeria monocytogenes", "Staphylococcus aureus",
    "Streptococcus pneumoniae", "Escherichia coli", "Campylobacter coli", "Campylobacter jejuni",
    "Legionella pneumophila", "Haemophilus influenzae", "Corynebacterium diphtheriae",
    "Vibrio alginolyticus", "Vibrio cholerae", "Bordetella pertussis", "Neisseria meningitidis", "Neisseria gonorrhoeae",
    "Shigella flexneri", "Shigella sonnei", "Shigella boydii", "Shigella dysenteriae",
]

# Exact rewrites of names in the sending file (lower case key). Rarely needed.
SPECIES_ALIASES = {
    "e.coli": "Escherichia coli",
}

# Detected values that mean "no species called" (case-insensitive). They count as a mismatch.
NO_CALL = {"", "unknown", "na", "n/a", "nan", "none", "??", "?"}

# Known confusions of the scanner: (expected genus, detected species), lower case.
# They are not counted as errors: the sample is shown as "unsure" (blue) instead of red.
CONFUSED_PAIRS = {("shigella", "escherichia coli")}

# Fill colours (ARGB)
COLOR_OK, COLOR_GENUS, COLOR_BAD, COLOR_NORUN = "FFA9D08E", "FFFFC000", "FFFF6B6B", "FFD9D9D9"
COLOR_BAD_LIGHT, COLOR_GENUS_LIGHT, COLOR_HEADER = "FFFFC7CE", "FFFFE699", "FFE7E6E6"
COLOR_UNSURE, COLOR_UNSURE_LIGHT = "FF9DC3E6", "FFDDEBF7"

# ====================================================================================
# Small helpers
# ====================================================================================
WELL_RE = re.compile(r"^\s*([A-Ha-h])\s*0*([1-9]|1[0-2])\s*$")
_XML_BAD = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
M = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def collapse(s):
    return re.sub(r"\s+", " ", str(s)).strip()


def norm_header(s):
    s = unicodedata.normalize("NFKD", str(s))
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    return re.sub(r"\s+", " ", s.lower()).strip().rstrip(".:# ")


def cell_text(v):
    if v is None:
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return collapse(v)


def parse_plate(v):
    if v is None or (isinstance(v, str) and not v.strip()):
        return None
    if isinstance(v, float) and v.is_integer():
        return int(v)
    if isinstance(v, str):
        s = v.strip()
        return int(s) if re.fullmatch(r"\d+(\.0+)?", s) else s
    return v


def plate_key(p):
    return (0, p, "") if isinstance(p, int) else (1, 0, str(p))


def col_letter(n):
    s = ""
    while n:
        n, rem = divmod(n - 1, 26)
        s = chr(65 + rem) + s
    return s


# ====================================================================================
# xlsx reader (stdlib)
# ====================================================================================
def _col_index(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n


def read_xlsx(path):
    """Return {sheet name: rows}, rows = list of lists of values (index 0 = column A)."""
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        shared = []
        if "xl/sharedStrings.xml" in names:
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("{%s}si" % M):
                text = "".join(t.text or "" for t in si.findall("{%s}t" % M))
                for r in si.findall("{%s}r" % M):          # rich text runs
                    t = r.find("{%s}t" % M)
                    text += (t.text or "") if t is not None else ""
                shared.append(text)

        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = {}
        if "xl/_rels/workbook.xml.rels" in names:
            for rel in ET.fromstring(z.read("xl/_rels/workbook.xml.rels")):
                rels[rel.get("Id")] = rel.get("Target")

        book = {}
        for sh in wb.find("{%s}sheets" % M):
            target = rels.get(sh.get("{%s}id" % R))
            if not target:
                continue
            target = target[1:] if target.startswith("/") else "xl/" + target
            if target not in names:
                continue
            grid, last_r = {}, 0
            for row in ET.fromstring(z.read(target)).iter("{%s}row" % M):
                r = int(row.get("r") or last_r + 1)
                last_r, last_c = r, 0
                for c in row.findall("{%s}c" % M):
                    ref = c.get("r")
                    ci = _col_index(ref) if ref else last_c + 1
                    last_c = ci
                    t, v, is_ = c.get("t"), c.find("{%s}v" % M), c.find("{%s}is" % M)
                    val = None
                    if t == "inlineStr" and is_ is not None:
                        val = "".join(x.text or "" for x in is_.iter("{%s}t" % M))
                    elif v is None or v.text is None:
                        val = None
                    elif t == "s":
                        val = shared[int(v.text)]
                    elif t in ("str", "d"):
                        val = v.text
                    elif t == "b":
                        val = v.text == "1"
                    elif t == "e":
                        val = None
                    else:
                        try:
                            f = float(v.text)
                            val = int(f) if f.is_integer() else f
                        except ValueError:
                            val = v.text
                    if val is not None:
                        grid[(r, ci)] = val
            n_rows = max((r for r, _ in grid), default=0)
            n_cols = max((c for _, c in grid), default=0)
            book[sh.get("name")] = [[grid.get((r, c)) for c in range(1, n_cols + 1)]
                                    for r in range(1, n_rows + 1)]
        return book


# ====================================================================================
# Sending file
# ====================================================================================
def find_table(rows):
    """The table can start on any row. Return (header row index, {field: column index})."""
    for i, row in enumerate(rows):
        cols = {}
        for j, cell in enumerate(row):
            if cell is None:
                continue
            h = norm_header(cell)
            for field, aliases in HEADER_ALIASES.items():
                if field not in cols and h in aliases:
                    cols[field] = j
        if len(cols) == len(HEADER_ALIASES) and len(set(cols.values())) == len(cols):
            return i, cols
    return None


def parse_sending(path, sheet=None):
    book = read_xlsx(path)
    if sheet and sheet not in book:
        raise SystemExit(f"Sheet {sheet!r} not found. Sheets: {list(book)}")
    for name in ([sheet] if sheet else list(book)):
        found = find_table(book[name])
        if found:
            break
    else:
        wanted = {k: v[0] for k, v in HEADER_ALIASES.items()}
        raise SystemExit("Could not find the table in the sending file: no row contains the columns "
                         f"{wanted} (see HEADER_ALIASES at the top of the script).")
    rows, (hdr, cols) = book[name], found
    samples, issues = [], []
    for i in range(hdr + 1, len(rows)):
        row = rows[i]

        def get(field):
            j = cols[field]
            return row[j] if j < len(row) else None

        code, sp = cell_text(get("code")), cell_text(get("species"))
        plate, well = parse_plate(get("plate")), cell_text(get("well"))
        if not code:
            if sp:
                issues.append(f"Row {i + 1}: species '{sp}' without a sample code (plate {plate}, well {well})")
            continue  # empty well
        samples.append({"code": code, "species": sp, "plate": plate, "well": well, "row": i + 1})
    return name, hdr + 1, samples, issues


# ====================================================================================
# Results (Terra export)
# ====================================================================================
def peek_columns(path):
    """Header of a results file (TSV/CSV), e.g. to let a UI offer the columns."""
    with open(path, encoding="utf-8-sig", newline="") as f:
        first = f.readline()
    delim = "\t" if "\t" in first else (";" if first.count(";") > first.count(",") else ",")
    return [h.strip() for h in next(csv.reader([first], delimiter=delim), [])]


def read_results(paths, id_col=None, species_col=None):
    """Accepted formats (columns auto-detected, or forced with id_col / species_col):
       - Terra export:  entity:xxx_id, ..., sp_detection_species_detected
       - Report:        Sample,Detected,Match   (the Match column is ignored: it is recomputed here)"""
    res, issues = {}, []
    for path in paths:
        with open(path, encoding="utf-8-sig", newline="") as f:
            text = f.read()
        first = text.splitlines()[0] if text.strip() else ""
        delim = "\t" if "\t" in first else (";" if first.count(";") > first.count(",") else ",")
        rows = list(csv.reader(io.StringIO(text), delimiter=delim))
        if not rows:
            continue
        header = [h.strip() for h in rows[0]]
        low = [h.lower() for h in header]

        def pick(name, candidates):
            if name:
                if name not in header:
                    raise ValueError(f"{path}: column {name!r} not found. Columns: {header}")
                return header.index(name)
            for c in candidates:
                for j, h in enumerate(low):
                    if c(h):
                        return j
            return None

        ci = pick(id_col, [lambda h: h.startswith("entity:"),
                           lambda h: h in ("sample_id", "sample", "id", "sample_name")])
        ci = 0 if ci is None else ci
        si = pick(species_col, [lambda h: h.endswith("species_detected"),
                                lambda h: h in ("detected_species", "detected", "species")])
        if si is None:
            raise ValueError(f"{path}: could not find the detected-species column (looked for "
                             f"'...species_detected', 'detected_species', 'detected'). Columns: {header}. "
                             f"Use --species-col.")
        for r in rows[1:]:
            if not r or ci >= len(r) or not r[ci].strip():
                continue
            sid = r[ci].strip()
            if sid in res:
                issues.append(f"Result for '{sid}' appears several times: the last one is used")
            res[sid] = collapse(r[si]) if si < len(r) else ""
    return res, issues


# ====================================================================================
# Species names
# ====================================================================================
class Namer:
    """Normalises names written by hand in the sending file ('S. aureus', 'e.coli', typos)."""

    def __init__(self, extra_names=()):
        self.pool = {}
        for sp in list(KNOWN_SPECIES) + [n for n in extra_names if n.lower() not in NO_CALL]:
            self.pool.setdefault(sp.lower(), sp)
            self.pool.setdefault(sp.split()[0].lower(), sp.split()[0])
        self.genera = sorted({sp.split()[0] for sp in KNOWN_SPECIES})
        self.aliases = {k.lower(): v for k, v in SPECIES_ALIASES.items()}
        self.changes, self.unknown = {}, set()

    def _expand(self, low):
        m = re.match(r"^([a-z])\.\s*([a-z][a-z\-]*)", low)  # "s. aureus", "e.coli", "v. chol non-o1"
        if not m:
            return None
        initial, tok = m.groups()
        cands = {sp for sp in KNOWN_SPECIES
                 if sp[0].lower() == initial and sp.split()[1].lower().startswith(tok)}
        if len(cands) == 1:
            return cands.pop()
        gen = [g for g in self.genera if g[0].lower() == initial]
        return f"{gen[0]} {tok}" if len(gen) == 1 else None

    def canon(self, raw):
        s = collapse(raw)
        if not s:
            return ""
        low = s.lower()
        out = self.aliases.get(low) or self.pool.get(low) or self._expand(low)
        if out is None:
            close = difflib.get_close_matches(low, list(self.pool), n=1, cutoff=0.9)
            out = self.pool[close[0]] if close else s
        if out.lower() not in self.pool:
            self.unknown.add(out)
        if out != s:
            self.changes[s] = out
        return out


def judge(expected, detected):
    """'ok' = match (same species, or same genus when only a genus is expected),
    'genus' = same genus but different species, 'unsure' = a known confusion of the scanner
    (CONFUSED_PAIRS), 'bad' = anything else."""
    e, d = expected.lower().split(), detected.lower().split()
    if detected.lower() in NO_CALL or not d:
        return "bad"
    if len(e) == 1:
        res = "ok" if d[0] == e[0] else "bad"
    elif e == d:
        res = "ok"
    else:
        res = "genus" if d[0] == e[0] else "bad"
    if res == "bad" and (e[0], " ".join(d)) in CONFUSED_PAIRS:
        return "unsure"
    return res


# ====================================================================================
# xlsx writer (stdlib)
# ====================================================================================
# style ids (index in cellXfs below)
S_DEF, S_TITLE, S_HEAD, S_OK, S_GENUS, S_BAD, S_NORUN, S_EMPTY = range(8)
S_ERRROW, S_ERRHEAD, S_BOLD, S_WRAP, S_GENUSROW, S_GENUSHEAD = range(8, 14)
S_UNSURE, S_UNSUREROW, S_UNSUREHEAD = range(14, 17)


def styles_xml():
    # explicit black: with "automatic" colour some viewers switch to white text on the red fill
    fonts = ['<font><sz val="11"/><color rgb="FF000000"/><name val="Calibri"/></font>',
             '<font><b/><sz val="11"/><color rgb="FF000000"/><name val="Calibri"/></font>',
             '<font><b/><sz val="14"/><color rgb="FF000000"/><name val="Calibri"/></font>',
             '<font><i/><sz val="11"/><color rgb="FF595959"/><name val="Calibri"/></font>',
             '<font><sz val="11"/><color rgb="FF9C0006"/><name val="Calibri"/></font>',
             '<font><b/><sz val="11"/><color rgb="FF9C0006"/><name val="Calibri"/></font>']

    def fill(rgb):
        return f'<fill><patternFill patternType="solid"><fgColor rgb="{rgb}"/><bgColor indexed="64"/></patternFill></fill>'

    fills = ['<fill><patternFill patternType="none"/></fill>', '<fill><patternFill patternType="gray125"/></fill>',
             fill(COLOR_HEADER), fill(COLOR_OK), fill(COLOR_GENUS), fill(COLOR_BAD), fill(COLOR_NORUN),
             fill(COLOR_BAD_LIGHT), fill(COLOR_GENUS_LIGHT), fill(COLOR_UNSURE), fill(COLOR_UNSURE_LIGHT)]
    side = '<{0} style="thin"><color rgb="FFBFBFBF"/></{0}>'
    borders = ['<border><left/><right/><top/><bottom/><diagonal/></border>',
               "<border>" + "".join(side.format(x) for x in ("left", "right", "top", "bottom")) + "<diagonal/></border>"]
    center = '<alignment horizontal="center" vertical="center"/>'
    left = '<alignment horizontal="left" vertical="center"/>'
    wrap = '<alignment vertical="top" wrapText="1"/>'
    # (font, fill, border, alignment)
    xfs = [(0, 0, 0, None), (2, 0, 0, None), (1, 2, 1, center), (0, 3, 1, center), (0, 4, 1, center),
           (0, 5, 1, center), (3, 6, 1, center), (0, 0, 1, center),
           (4, 7, 1, left), (5, 7, 1, left), (1, 0, 0, None), (0, 0, 0, wrap), (0, 8, 1, left), (1, 8, 1, left),
           (0, 9, 1, center), (0, 10, 1, left), (1, 10, 1, left)]
    xf_xml = []
    for fo, fi, bo, al in xfs:
        attrs = f'numFmtId="0" fontId="{fo}" fillId="{fi}" borderId="{bo}" xfId="0" applyFont="1" applyFill="1" applyBorder="1"'
        xf_xml.append(f"<xf {attrs} applyAlignment=\"1\">{al}</xf>" if al else f"<xf {attrs}/>")
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            f'<styleSheet xmlns="{M}">'
            f'<fonts count="{len(fonts)}">{"".join(fonts)}</fonts>'
            f'<fills count="{len(fills)}">{"".join(fills)}</fills>'
            f'<borders count="{len(borders)}">{"".join(borders)}</borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            f'<cellXfs count="{len(xf_xml)}">{"".join(xf_xml)}</cellXfs>'
            '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
            '</styleSheet>')


class Sheet:
    def __init__(self, name, tab=None, grid=True):
        self.name, self.tab, self.grid = name, tab, grid
        self.cells, self.widths, self.heights, self.merges = {}, {}, {}, []

    def put(self, r, c, v=None, s=S_DEF):
        self.cells[(r, c)] = (v, s)

    def merge(self, r, c1, c2, v, s):
        self.put(r, c1, v, s)
        for c in range(c1 + 1, c2 + 1):
            self.put(r, c, None, s)
        self.merges.append(f"{col_letter(c1)}{r}:{col_letter(c2)}{r}")

    def xml(self, selected=False):
        rows = {}
        for (r, c), (v, s) in self.cells.items():
            rows.setdefault(r, []).append((c, v, s))
        out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', f'<worksheet xmlns="{M}">']
        tab = f'<tabColor rgb="{self.tab}"/>' if self.tab else ""
        out.append(f'<sheetPr>{tab}<pageSetUpPr fitToPage="1"/></sheetPr>')   # print on one page width
        out.append('<sheetViews><sheetView workbookViewId="0"%s%s/></sheetViews>'
                   % ("" if self.grid else ' showGridLines="0"', ' tabSelected="1"' if selected else ""))
        out.append('<sheetFormatPr defaultRowHeight="15"/>')
        if self.widths:
            out.append("<cols>" + "".join(
                f'<col min="{c}" max="{c}" width="{w}" customWidth="1"/>' for c, w in sorted(self.widths.items())) + "</cols>")
        out.append("<sheetData>")
        for r in sorted(rows):
            ht = f' ht="{self.heights[r]}" customHeight="1"' if r in self.heights else ""
            out.append(f'<row r="{r}"{ht}>')
            for c, v, s in sorted(rows[r], key=lambda x: x[0]):
                ref = f"{col_letter(c)}{r}"
                if v is None:
                    out.append(f'<c r="{ref}" s="{s}"/>')
                elif isinstance(v, (int, float)) and not isinstance(v, bool):
                    out.append(f'<c r="{ref}" s="{s}"><v>{v}</v></c>')
                else:
                    txt = escape(_XML_BAD.sub("", str(v)))
                    out.append(f'<c r="{ref}" s="{s}" t="inlineStr"><is><t xml:space="preserve">{txt}</t></is></c>')
            out.append("</row>")
        out.append("</sheetData>")
        if self.merges:
            out.append(f'<mergeCells count="{len(self.merges)}">'
                       + "".join(f'<mergeCell ref="{m}"/>' for m in self.merges) + "</mergeCells>")
        out.append('<pageMargins left="0.5" right="0.5" top="0.6" bottom="0.6" header="0.3" footer="0.3"/>')
        out.append('<pageSetup orientation="landscape" fitToWidth="1" fitToHeight="0"/>')
        out.append("</worksheet>")
        return "".join(out)


def write_xlsx(path, sheets):
    names, used = [], set()
    for sh in sheets:                                   # valid + unique sheet names (max 31 chars)
        base = re.sub(r"[\[\]:*?/\\]", "-", sh.name)[:31] or "Sheet"
        n, k = base, 2
        while n.lower() in used:
            n = f"{base[:28]}-{k}"
            k += 1
        used.add(n.lower())
        names.append(n)
    ct = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
          '<Default Extension="xml" ContentType="application/xml"/>'
          '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
          '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
          + "".join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                    for i in range(1, len(sheets) + 1)) + "</Types>")
    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 f'<Relationship Id="rId1" Type="{R}/officeDocument" Target="xl/workbook.xml"/></Relationships>')
    wb = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          f'<workbook xmlns="{M}" xmlns:r="{R}"><bookViews><workbookView/></bookViews><sheets>'
          + "".join(f'<sheet name="{escape(n, {chr(34): "&quot;"})}" sheetId="{i}" r:id="rId{i}"/>'
                    for i, n in enumerate(names, 1)) + "</sheets></workbook>")
    wb_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
               '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
               + "".join(f'<Relationship Id="rId{i}" Type="{R}/worksheet" Target="worksheets/sheet{i}.xml"/>'
                         for i in range(1, len(sheets) + 1))
               + f'<Relationship Id="rId{len(sheets) + 1}" Type="{R}/styles" Target="styles.xml"/></Relationships>')
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", ct)
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        z.writestr("xl/styles.xml", styles_xml())
        for i, sh in enumerate(sheets, 1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", sh.xml(selected=(i == 1)))


# ====================================================================================
# Build the report
# ====================================================================================
STYLE_OF = {"ok": S_OK, "genus": S_GENUS, "bad": S_BAD, "norun": S_NORUN, "unsure": S_UNSURE}


def plate_sheet(plate, items, counts):
    tab = (COLOR_BAD if counts["bad"] else COLOR_GENUS if counts["genus"]
           else COLOR_UNSURE if counts["unsure"] else COLOR_OK)
    sh = Sheet(f"Plate {plate}", tab=tab, grid=False)
    sh.widths = {c: 14 for c in range(1, 14)}
    sh.put(1, 1, f"Plate No. {plate}", S_TITLE)
    sh.put(2, 1, "Legend:", S_BOLD)
    for c, (txt, st) in enumerate([("match", S_OK), ("same genus", S_GENUS), ("mismatch", S_BAD),
                                   ("unsure", S_UNSURE), ("no result", S_NORUN)], start=2):
        sh.put(2, c, txt, st)
    sh.put(4, 1, None, S_HEAD)
    for c in range(1, 13):
        sh.put(4, 1 + c, c, S_HEAD)
    for r, letter in enumerate("ABCDEFGH"):
        sh.put(5 + r, 1, letter, S_HEAD)
        sh.heights[5 + r] = 24
        for c in range(1, 13):
            sh.put(5 + r, 1 + c, None, S_EMPTY)
    taken = {}
    for it in items:
        m = WELL_RE.match(it["well"])
        if m and (m.group(1).upper(), int(m.group(2))) not in taken:
            r, c = "ABCDEFGH".index(m.group(1).upper()), int(m.group(2))
            taken[(m.group(1).upper(), c)] = it
            sh.put(5 + r, 1 + c, it["code"], STYLE_OF[it["status"]])

    row = 14

    def listing(title, status, head_style, row_style, extra_note=""):
        nonlocal row
        sel = [it for it in items if it["status"] == status]
        if not sel:
            return
        sh.merge(row, 1, 12, f"{title} ({len(sel)}){extra_note}", head_style)
        row += 1
        sh.put(row, 1, "Sample ID", head_style)
        sh.merge(row, 2, 5, "Detected", head_style)
        sh.merge(row, 6, 9, "Expected (sending file)", head_style)
        sh.put(row, 10, "Well", head_style)
        row += 1
        for it in sorted(sel, key=lambda x: (x["well_sort"], x["code"])):
            sh.put(row, 1, it["code"], row_style)
            sh.merge(row, 2, 5, it["detected"] or "(none)", row_style)
            sh.merge(row, 6, 9, it["expected"] or "(none)", row_style)
            sh.put(row, 10, it["well"], row_style)
            row += 1
        row += 1

    listing("Errors: detected species does not match", "bad", S_ERRHEAD, S_ERRROW)
    listing("Same genus, different species", "genus", S_GENUSHEAD, S_GENUSROW)
    listing("Unsure: known confusion of the scanner (e.g. Shigella detected as E. coli)", "unsure",
            S_UNSUREHEAD, S_UNSUREROW)
    return sh


def build_report(sending, results, hdr_info, send_issues, res_issues, out_path):
    sheet_name, hdr_row, samples = hdr_info
    namer = Namer(extra_names=set(results.values()))
    res_lower = {k.lower(): (k, v) for k, v in results.items()}
    notes, used_ids = [], set()

    plates = {}
    for s in samples:
        s["expected"] = namer.canon(s["species"])
        hit = res_lower.get(s["code"].lower())
        if hit is None and len(s["code"]) >= 4:        # result id = code + suffix (e.g. SH037528_S12)
            k = s["code"].lower()
            cands = [v for kk, v in res_lower.items()
                     if kk.startswith(k) and len(kk) > len(k) and not kk[len(k)].isalnum()]
            if len(cands) == 1:
                hit = cands[0]
                notes.append(f"Result id '{hit[0]}' was matched to sample code '{s['code']}' (prefix match)")
        if hit is None:
            s["status"], s["detected"] = "norun", None
        else:
            used_ids.add(hit[0])
            s["detected"] = hit[1]
            s["status"] = judge(s["expected"], hit[1]) if s["expected"] else "norun"
            if not s["expected"]:
                notes.append(f"{s['code']}: no species in the sending file, cannot be checked (grey)")
        m = WELL_RE.match(s["well"])
        s["well_sort"] = ("ABCDEFGH".index(m.group(1).upper()) + 8 * int(m.group(2))) if m else 999
        plates.setdefault(s["plate"], []).append(s)

    sheets, summary, any_result = [], [], False
    for plate in sorted(plates, key=plate_key):
        items = plates[plate]
        counts = {k: sum(1 for i in items if i["status"] == k) for k in ("ok", "genus", "bad", "unsure", "norun")}
        n_res = len(items) - counts["norun"]
        if plate is None:
            for i in items:
                if i["detected"] is not None:
                    notes.append(f"{i['code']} has a result but no plate in the sending file")
            continue
        if n_res == 0:
            continue                                    # no sample of this plate was run: no sheet
        any_result = True
        seen = {}
        for i in items:
            m = WELL_RE.match(i["well"])
            if not m:
                notes.append(f"Plate {plate}: sample {i['code']} has an invalid well '{i['well']}' (not drawn in the grid)")
            else:
                key = (m.group(1).upper(), int(m.group(2)))
                if key in seen:
                    notes.append(f"Plate {plate}: well {i['well']} is used by {seen[key]} and {i['code']} (only the first is drawn)")
                seen.setdefault(key, i["code"])
        sheets.append(plate_sheet(plate, items, counts))
        summary.append((plate, len(items), n_res, counts))

    if not any_result:
        ex_s = ", ".join(s["code"] for s in samples[:3])
        ex_r = ", ".join(list(results)[:3])
        raise SystemExit("None of the result ids matches a sample code of the sending file.\n"
                         f"  sending file (sheet '{sheet_name}', header row {hdr_row}): {ex_s} ...\n"
                         f"  results: {ex_r} ...\n  Check --id-col.")

    n = Sheet("Notes")
    n.widths = {1: 120}
    r = 1
    n.put(r, 1, "Summary", S_TITLE)
    r += 1
    n.put(r, 1, f"Sending file: sheet '{sheet_name}', table header on row {hdr_row}, {len(samples)} samples", S_WRAP)
    r += 1
    for plate, tot, nres, c in summary:
        n.put(r, 1, f"Plate {plate}: {nres}/{tot} samples with a result: {c['ok']} match, "
                    f"{c['genus']} same genus, {c['bad']} mismatch, {c['unsure']} unsure", S_WRAP)
        r += 1
    sections = [
        ("Species names normalised from the sending file (check that they are right)",
         [f"{a}  ->  {b}" for a, b in sorted(namer.changes.items())]),
        ("Species names not in the known list (compared as typed; add them to KNOWN_SPECIES if needed)",
         sorted(namer.unknown)),
        ("Results that are not in the sending file",
         [f"{k}: {v or '(none)'}" for k, v in results.items() if k not in used_ids]),
        ("Warnings", send_issues + res_issues + notes),
    ]
    for title, lines in sections:
        if not lines:
            continue
        r += 1
        n.put(r, 1, title, S_BOLD)
        for ln in lines:
            r += 1
            n.put(r, 1, ln, S_WRAP)
    sheets.append(n)
    write_xlsx(out_path, sheets)
    return summary, namer, [l for _, ls in sections for l in ls]


# ====================================================================================
# CLI
# ====================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description="Check sequencing results against the sending file, plate by plate.")
    ap.add_argument("--sending", required=True, help="Sending file (.xlsx)")
    ap.add_argument("--results", required=True, nargs="+", help="Results TSV/CSV: Terra export, or a report Sample,Detected,Match")
    ap.add_argument("--out", default="plate_check.xlsx", help="Output xlsx (default: plate_check.xlsx)")
    ap.add_argument("--sheet", help="Sheet of the sending file (default: first sheet with the table)")
    ap.add_argument("--id-col", help="Results column with the sample id (default: auto)")
    ap.add_argument("--species-col", help="Results column with the detected species (default: auto-detected)")
    args = ap.parse_args(argv)

    sheet_name, hdr_row, samples, send_issues = parse_sending(args.sending, args.sheet)
    try:
        results, res_issues = read_results(args.results, args.id_col, args.species_col)
    except ValueError as e:
        raise SystemExit(str(e))
    print(f"Sending file: sheet '{sheet_name}', header on row {hdr_row}, {len(samples)} samples")
    print(f"Results: {len(results)} samples")
    summary, namer, _ = build_report(args.sending, results, (sheet_name, hdr_row, samples),
                                     send_issues, res_issues, args.out)
    for plate, tot, nres, c in summary:
        print(f"  Plate {plate}: {nres}/{tot} with result | match {c['ok']} | same genus {c['genus']} | mismatch {c['bad']} | unsure {c['unsure']}")
    for a, b in sorted(namer.changes.items()):
        print(f"  name normalised: '{a}' -> '{b}'")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
