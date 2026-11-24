#!/usr/bin/env python3
"""
Patch meningotype v0.8.6b pour Python 3.

Usage: python patch_meningotype.py /chemin/vers/meningotype_src/meningotype
       (le dossier du package, celui qui contient meningotype.py)

Le script échoue (assert) si un correctif ne s'applique pas,
pour que le build Docker s'arrête au lieu de produire une image cassée.
"""
import pathlib
import re
import sys

if len(sys.argv) != 2:
    sys.exit("Usage: patch_meningotype.py <dossier_du_package>")

pkg = pathlib.Path(sys.argv[1])
for name in ("meningotype.py", "ctrA.py", "run_blast.py"):
    assert (pkg / name).exists(), "fichier introuvable: %s" % (pkg / name)

# ---------------------------------------------------------------------------
# 1. run_blast.py : remplacement complet
#    (accepte seg/matrix/query_gencode/ungapped/..., nettoie les outfmt)
# ---------------------------------------------------------------------------
RUN_BLAST = r'''import subprocess
import sys

# Options BLAST qui sont des switches (sans valeur)
_SWITCHES = {"ungapped"}


def _clean(value):
    """Retire les guillemets/apostrophes litteraux (code ecrit pour shell=True)."""
    return str(value).strip("'\"")


def seqBLAST(query, db, blast, outfmt, perc_identity=None, evalue=None,
             num_threads=1, culling_limit=None, fasta_data=None, **kwargs):
    command = [blast, "-query", query, "-db", db, "-outfmt", _clean(outfmt)]

    if evalue is not None:
        command += ["-evalue", str(evalue)]
    # -perc_identity n'existe que pour blastn
    if perc_identity is not None and blast == "blastn":
        command += ["-perc_identity", str(perc_identity)]
    if culling_limit is not None:
        command += ["-culling_limit", str(culling_limit)]
    command += ["-num_threads", str(num_threads)]

    # Options supplementaires (task, seg, matrix, query_gencode, comp_based_stats, ungapped...)
    for key, value in kwargs.items():
        if value is None:
            continue
        if key in _SWITCHES:
            if str(value).lower() in ("true", "yes", "1"):
                command.append("-" + key)
        else:
            command += ["-" + key, _clean(value)]

    try:
        result = subprocess.run(
            command,
            input=fasta_data if query == "-" else None,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        sys.stderr.write("BLAST failed: %s\n%s\n" % (" ".join(command), e.stderr))
        raise
    return result.stdout, result.stderr
'''
(pkg / "run_blast.py").write_text(RUN_BLAST)

# ---------------------------------------------------------------------------
# 2. meningotype.py : urllib, https, bytes -> str
# ---------------------------------------------------------------------------
p = pkg / "meningotype.py"
s = p.read_text()

# urllib.urlretrieve n'existe plus en py3
s = s.replace("urllib.urlretrieve", "urllib.request.urlretrieve")

# imports (apres __future__, qui doit rester en premier)
if "import urllib.request" not in s:
    s = s.replace(
        "from __future__ import print_function",
        "from __future__ import print_function\n"
        "import urllib.request\n"
        "from urllib.error import HTTPError",
        1,
    )

# API PubMLST en https
s = s.replace("http://rest.pubmlst.org", "https://rest.pubmlst.org")

# Popen.communicate() renvoie des bytes en py3 (mlst / isPcr)
s = re.sub(
    r"(proc\.communicate\(\)\[0\])[ \t]*$",
    r"\1.decode('utf-8')",
    s,
    flags=re.M,
)

p.write_text(s)

assert "urllib.urlretrieve" not in s, "patch urllib: urlretrieve encore present"
assert "import urllib.request" in s, "patch urllib: import manquant"
assert "http://rest.pubmlst.org" not in s, "patch https non applique"
assert re.search(r"proc\.communicate\(\)\[0\][ \t]*$", s, flags=re.M) is None, \
    "patch decode: un communicate() n'est pas decode"

# ---------------------------------------------------------------------------
# 3. ctrA.py : NcbiblastnCommandline n'est plus importe -> run_blast.seqBLAST
# ---------------------------------------------------------------------------
c = pkg / "ctrA.py"
cs = c.read_text()

cs, n = re.subn(
    r"ctrABLAST = NcbiblastnCommandline\(.*?stdout, stderr = ctrABLAST\(\)",
    "stdout, stderr = run_blast.seqBLAST(query=f, db=ctrADB, blast='blastn', "
    "task='blastn', outfmt='6 sseqid pident length', perc_identity=90, "
    "evalue='1e-20', culling_limit=1)",
    cs,
    flags=re.S,
)
assert n == 1, "patch ctrA.py: bloc NcbiblastnCommandline introuvable"
c.write_text(cs)

assert not re.search(r"^[^#\n]*NcbiblastnCommandline\(", cs, flags=re.M), \
    "ctrA.py: appel NcbiblastnCommandline encore present"

# ---------------------------------------------------------------------------
# 4. Verification de syntaxe
# ---------------------------------------------------------------------------
import py_compile
for f in sorted(pkg.glob("*.py")):
    py_compile.compile(str(f), doraise=True)

print("meningotype patched OK")