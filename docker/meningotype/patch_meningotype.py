import pathlib
import re
import sys

pkg = pathlib.Path(sys.argv[1])  # dossier du package (celui qui contient meningotype.py)

RUN_BLAST = '''import subprocess
import sys

# Options BLAST qui sont des switches (sans valeur)
_SWITCHES = {"ungapped"}


def _clean(value):
    """Retire les guillemets/apostrophes litteraux (code ecrit pour shell=True)."""
    return str(value).strip("'\\"")


def seqBLAST(query, db, blast, outfmt, perc_identity=None, evalue=None,
             num_threads=1, culling_limit=None, fasta_data=None, **kwargs):
    command = [blast, "-query", query, "-db", db, "-outfmt", _clean(outfmt)]

    if evalue is not None:
        command += ["-evalue", str(evalue)]
    if perc_identity is not None and blast == "blastn":
        command += ["-perc_identity", str(perc_identity)]
    if culling_limit is not None:
        command += ["-culling_limit", str(culling_limit)]
    command += ["-num_threads", str(num_threads)]

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
        sys.stderr.write("BLAST failed: %s\\n%s\\n" % (" ".join(command), e.stderr))
        raise
    return result.stdout, result.stderr
'''

# 1. run_blast.py
(pkg / "run_blast.py").write_text(RUN_BLAST)

# 2. meningotype.py
p = pkg / "meningotype.py"
s = p.read_text()

s = s.replace("urllib.urlretrieve", "urllib.request.urlretrieve")
if "import urllib.request" not in s:
    s = s.replace(
        "from __future__ import print_function",
        "from __future__ import print_function\nimport urllib.request\nfrom urllib.error import HTTPError",
        1,
    )
s = s.replace("http://rest.pubmlst.org", "https://rest.pubmlst.org")
# Popen.communicate() renvoie des bytes en py3 (mlst / isPcr)
s = re.sub(r"(proc\.communicate\(\)\[0\])[ \t]*$", r"\1.decode('utf-8')", s, flags=re.M)

p.write_text(s)

# Verifications : echec du build si un patch n'a pas pris
assert "urllib.request.urlretrieve" in s, "patch urllib manquant"
assert "import urllib.request" in s, "import urllib.request manquant"
assert "communicate()[0].decode('utf-8')" in s, "patch decode manquant"
print("meningotype patched OK")