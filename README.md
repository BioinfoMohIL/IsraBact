# 🧬 IsraBact WDL Workflows

<p align="center">
  <img src="https://img.shields.io/badge/WDL-1.0-4B8BBE?style=for-the-badge&logo=python" alt="WDL 1.0">
  <img src="https://img.shields.io/badge/Docker-Ready-2496ED?style=for-the-badge&logo=docker" alt="Docker">
  <img src="https://img.shields.io/badge/Status-Production-00C853?style=for-the-badge" alt="Status">
  <img src="https://img.shields.io/badge/Israel-MOH-0038B8?style=for-the-badge" alt="Israel MOH">
</p>

<p align="center">
  <b>🦠 Complete bacterial genomics toolkit for public health surveillance</b><br>
  From raw reads to virulence prediction in one pipeline
</p>

---

## 🎯 What's Inside

| Workflow | What It Does | Input | Output |
|----------|--------------|-------|--------|
| 🔬 **Species Detection** | Check if your samples are what you think they are | Basespace reads | Species match ✅/❌ |
| ☠️ **Virulence Finder** | Find nasty virulence genes | Assembly FASTA | Gene presence matrix |
| 🧫 **ListPred** | Predict *Listeria* serotype & phenotype | FASTA or FASTQ | Virulence class + disinfectant tolerance |
| 🧬 **fq2DNA Assembly** | Build genomes from scratch | Paired-end reads | Quality assembly + metrics |
| 🔧 **Long Reads Assembly** | Unicycler + Pilon polishing | Short + long reads | Hybrid assembly + QC |
| 🦠 **DiphtOscan** | *C. diphtheriae* comprehensive analysis | Assembly FASTA | MLST + toxin + AMR + tree |

---

## 🚀 Quick Start

### 1️⃣ Run Species Check

```json
{
  "api_server": "https://api.basespace.illumina.com",
  "access_token": "YOUR_TOKEN",
  "basespace_collection_id": "12345",
  "sample_prefix": "EC"
}
```

### 2️⃣ Find Virulence Genes

```json
{
  "fasta": "genome.fasta",
  "species": "salmonella",
  "samplename": "outbreak_001"
}
```

### 3️⃣ Assemble Genome

```json
{
  "read1": "sample_R1.fastq.gz",
  "read2": "sample_R2.fastq.gz",
  "samplename": "sample_001"
}
```

---

## 📦 Workflows Deep Dive

### 🔬 Species Detection
**File:** `metagenomics/wf_species_detection_bs_reads.wdl`

```
Basespace → Kraken2 → Species Match Report
```

**Perfect for:** QC before starting expensive analyses

---

### ☠️ Virulence Finder
**File:** `virulence/wf_virulence_finder.wdl`

```
Assembly → BLAST vs FDA DB → Gene Matrix (CSV + HTML)
```

**Supports:** *Salmonella*, *E. coli*  
**Bonus:** Plasmid virulence check with `plasmid_check: true`

---

### 🧫 ListPred (Listeria)
**Files:** `resistance/wf_listpred_fasta.wdl` | `wf_listpred_reads.wdl`

```
FASTA/FASTQ → Snakemake → Serotype + Virulence + Disinfectant Tolerance
```

**Use case:** Food safety surveillance

---

### 🧬 fq2DNA Assembly
**File:** `utilities/wf_fq2dna.wdl`

```
Raw Reads → fq2DNA → Assembly + Metrics + Selected Scaffolds
```

**Special:** Auto-detects contamination with `alien_tag`

**Organism types:**
- `B` = Bacteria
- `V` = Virus
- `E` = Eukaryote

---

### 🔧 Long Reads Assembly
**File:** `utilities/wf_assembly_long_reads.wdl`

```
Short + Long Reads → Unicycler → Pilon Polishing → QC (QUAST + BUSCO)
```

**Perfect for:** Hybrid assemblies (Illumina + Nanopore/PacBio)

**Quality reports:**
- ✅ N50, assembly length
- ✅ GC content
- ✅ BUSCO completeness

---

### 🦠 DiphtOscan
**File:** `epidemiology/diphtheria/wf_diphtoscan.wdl`

```
Assembly → diphtOscan → MLST + Toxin + AMR + Tree
```

**Features:**
- 🧬 MLST typing
- 🦠 Toxin detection
- 💊 Resistance profiling
- 🌳 Phylogenetic tree (optional)

**[📖 Full docs →](DIPHTOSCAN_README.md)**

---

## 🎨 Output Examples

### Virulence Matrix (CSV)
```csv
Sample,gene_A,gene_B,gene_C
strain_1,1,0,1
strain_2,1,1,0
```

### Species Detection
```csv
Sample,Detected_Species,Expected,Match
EC001,Escherichia coli,Escherichia,+
EC002,Salmonella enterica,Escherichia,-
```

### ListPred Result
```
Virulence Class: HV (Hypervirulent)
Disinfectant Tolerance: Tolerant
```

---

## 🛠️ Requirements

- 🐳 **Docker** (containers for all tools)
- ⚙️ **MiniWDL** or **Cromwell** (workflow engine)
- ☁️ **Cloud** or **HPC** (optional but recommended)

### Install MiniWDL

```bash
pip install miniwdl
```

### Run a Workflow

```bash
miniwdl run wf_virulence_finder.wdl -i inputs.json
```

---

## 🎯 Use Cases

| Scenario | Recommended Workflow |
|----------|---------------------|
| 🔍 **Unknown sample QC** | Species Detection |
| 🦠 **Outbreak investigation** | Virulence Finder + DiphtOscan |
| 🧫 **Food safety** | ListPred |
| 🧬 **De novo assembly** | fq2DNA or Long Reads |
| 📊 **Surveillance** | All of the above! |

---

## 🏗️ Architecture

```
IsraBact WDL/
├── 📁 metagenomics/
│   └── wf_species_detection_bs_reads.wdl
├── 📁 virulence/
│   └── wf_virulence_finder.wdl
├── 📁 resistance/
│   ├── wf_listpred_fasta.wdl
│   └── wf_listpred_reads.wdl
├── 📁 utilities/
│   ├── wf_fq2dna.wdl
│   └── wf_assembly_long_reads.wdl
└── 📁 epidemiology/diphtheria/
    └── wf_diphtoscan.wdl
```

---

## 💡 Pro Tips

### Speed Up Assemblies
```json
{
  "cpu": 32,
  "call_pilon": false
}
```

### Batch Processing
```json
{
  "assemblies": ["s1.fasta", "s2.fasta", "s3.fasta"]
}
```

### Skip DB Updates
```json
{
  "update_db": false
}
```

---

## 🐛 Troubleshooting

### Docker not found?
```bash
docker --version
# Install: https://docs.docker.com/get-docker/
```

### Out of memory?
```json
{
  "memory": "32 GB"
}
```

### Workflow fails?
```bash
# Check logs
ls -la _LAST/
cat _LAST/stderr.txt
```

---

## 📚 Resources

- 🔧 **Theiagen Workflows**: [GitHub](https://github.com/theiagen/public_health_bioinformatics)
- 🦠 **FDA VirulenceFinder**: [virulence.fda.gov](https://virulence.fda.gov/)
- 🧫 **ListPred**: [GitHub](https://github.com/genomicepidemiology/ListPred)
- 🧬 **fq2DNA**: Institut Pasteur (Alexis Criscuolo)
- 🔬 **Kraken2**: [GitHub](https://github.com/DerrickWood/kraken2)
- 🦠 **DiphtOscan**: [GitLab Pasteur](https://gitlab.pasteur.fr/BEBP/diphtoscan)

---

## 👥 Team

**Developed by:**  
🏥 **Ministry of Health** – Jerusalem, Israel

**Maintainer:**  
👨‍💻 **David Maimoun**  
📧 david.maimoun@moh.gov.il

**Based on:**  
Theiagen Public Health Bioinformatics pipelines

---

## 📄 License

Ministry of Health – Jerusalem, Israel

---

## 🙏 Acknowledgments

Special thanks to:
- 🦠 Theiagen Genomics team
- 🧬 Institut Pasteur bioinformatics
- 🔬 FDA CFSAN researchers
- 🌍 Public health genomics community

---

## 🆘 Support

- 🐛 **Issues**: Open an issue on the repository
- 📧 **Email**: david.maimoun@moh.gov.il
- 💬 **Questions**: Tag maintainer in issues

---

<p align="center">
  <b>🧬 Built with ❤️ for public health genomics 🦠</b><br>
  <sub>Protecting communities through bacterial surveillance</sub>
</p>

---

## 📊 Quick Reference

| Need | Workflow | Time* |
|------|----------|-------|
| Species check | Species Detection | ~5 min |
| Virulence genes | Virulence Finder | ~10 min |
| *Listeria* typing | ListPred | ~15 min |
| De novo assembly | fq2DNA | ~30 min |
| Hybrid assembly | Long Reads | ~1-2 hrs |
| *Diphtheria* analysis | DiphtOscan | ~20 min |

*\*Approximate times for typical samples*

---

**Ready to start? Pick a workflow above and run it! 🚀**