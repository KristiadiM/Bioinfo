# 🧬 FoldBot — AI-Powered Protein Structure & In Silico Mutagenesis Partner

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Telegram Bot API](https://img.shields.io/badge/Telegram-Bot%20API-blue?logo=telegram)](https://core.telegram.org/bots)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ESMFold](https://img.shields.io/badge/Structure-ESMFold%20%7C%20AlphaFold-brightgreen)](https://esmatlas.com)
[![PoPMuSiC](https://img.shields.io/badge/Stability-PoPMuSiC%20%7C%20FoldX-orange)](https://soft.dezyme.com/)
[![LLM: Gemini + Groq](https://img.shields.io/badge/LLM-Gemini%20%2B%20Groq%20Fallback-purple)](https://groq.com/)

**FoldBot** is a Telegram bot that brings protein structure prediction, *in silico* mutagenesis, thermodynamic stability screening ($\Delta\Delta G$), and physicochemical analysis directly to your chat. 

Powered by **Meta ESMFold**, **AlphaFold DB**, **BioPython**, **PoPMuSiC-modeled statistical potentials**, and a dual-LLM reasoning engine (**Google Gemini** with automated **Groq Qwen 2.5** fallback), FoldBot transforms complex computational biology pipelines into natural language conversations and mobile-friendly outputs.

---

## Why I Created This Bot?

In modern life sciences, understanding how point mutations impact protein stability and 3D conformation is critical for drug discovery, enzyme engineering, and interpreting disease-associated genetic variants.

However, computational structural biology traditionally poses steep barriers to entry:
1. **High Hardware Requirements:** Running state-of-the-art predictors like AlphaFold or ESMFold locally requires workstations with massive GPU VRAM.
2. **Complex Software Environments:** Installing physics-based stability tools (FoldX, Rosetta, ThermoMPNN) often involves complex Linux dependencies, compiler configurations, and proprietary academic licenses.
3. **Fragmented Workflows:** A researcher typically has to jump between multiple separate tools: one server to predict structures, another to align PDBs, another to calculate $\Delta\Delta G$, and manual scripts to check hydropathy or charge changes.

### The Vision
I created **FoldBot** to **democratize access to structural bioinformatics**:
* **Accessible Anywhere:** Wet-lab biologists, students, and bioengineers can pull out their smartphone in the lab or classroom and evaluate protein mutations in under 30 seconds.
* **Zero Local Installation Required:** All heavy lifting is offloaded to cloud APIs and built-in analytical heuristics. You don't need a GPU cluster or local command-line tools to start exploring.
* **Natural Language Interpretation:** Rather than just outputting raw coordinate tables and energy matrices, the bot uses generative AI to explain *why* a mutation is destabilizing in plain, biochemically grounded language.

---

## Key Features

- **Instant 3D Structure Prediction:** 
  - Predict full atomic coordinates from raw amino acid sequences via the free **Meta ESMFold API** or enterprise **NVIDIA NIM ESMFold**.
  - Download experimental and predicted structures instantly via the **AlphaFold DB** (EMBL-EBI) with zero queue time.
- ***In Silico* Mutagenesis & Structure Superposition:**
  - Introduce single or multiple point mutations (e.g., `V11R`, `K11R,E20Q`).
  - Aligns mutant and wild-type structures using the **Kabsch algorithm** to calculate per-residue $\text{C}\alpha$ RMSD ($\text{Å}$) and backbone shifts.
  - Automatically renders dual-track **pLDDT confidence curves** and displacement plots.
- **Thermodynamic Stability ($\Delta\Delta G$) Analysis:**
  - **Built-in PoPMuSiC Statistical Potentials:** Predicts free energy changes ($\text{kcal/mol}$) based on solvent accessibility, core hydrophobic transfer energy, cavity vs. steric volume changes, and low-dielectric desolvation penalties—with **zero local installation**.
  - **Local Physics Engines (Optional):** Supports local executables for **FoldX 5.0** and **ThermoMPNN** when installed.
  - Direct links to the official [PoPMuSiC Web Server](https://soft.dezyme.com/) for external cross-validation.
- **Physicochemical Profiling:**
  - Comprehensive sequence characterization using **BioPython ProtParam**: Molecular weight, isoelectric point ($\text{pI}$), net charge at $\text{pH } 7.0$, GRAVY hydropathy, instability index, and extinction coefficients.
- **Resilient Dual-LLM Conversational Engine:**
  - Seamless natural language interaction: chat with the bot normally without needing slash commands.
  - **Automatic Failover:** Uses Google Gemini as primary reasoning; if Google's free tier encounters capacity limits (HTTP 503), the bot automatically and silently redirects prompts to **Qwen 2.5 on Groq** in $\sim 30\text{ ms}$.
- **Automated Molecular Visualization:**
  - Generates publication-ready 3D cartoon ribbon renderings colored by pLDDT confidence (blue $\ge 90$, cyan $70\text{--}90$, yellow $50\text{--}70$, orange $< 50$).
  - Sends full `.pdb` files directly to Telegram for immediate inspection in **Mol\***, **PyMOL**, or **ChimeraX**.

---

## Scientific Methodology & Scoring Reference

FoldBot is built on transparent, peer-reviewed computational biology methods. Every metric reported by the bot is grounded in specific biophysical scales:

| Analysis | Metric / Score | Scientific Algorithm & Reference |
| :--- | :--- | :--- |
| **3D Coordinate Prediction** | Structure (`.pdb`) | **ESMFold** (Lin et al., *Science* 2023) / **AlphaFold 2** (Jumper et al., *Nature* 2021 via EMBL-EBI) |
| **Model Confidence** | `pLDDT` (0–100 scale) | **predicted Local Distance Difference Test** (Mariani et al., *Bioinformatics* 2013) |
| **Structural Alignment** | $\text{C}\alpha\text{ RMSD}$ & Backbone Shift ($\text{Å}$) | **Kabsch 3D coordinate superposition** (Kabsch, *Acta Cryst.* 1976 via BioPython `Superimposer`) |
| **Residue Burial / Environment** | Core vs. Surface | **$10\text{ Å}$ $\text{C}\alpha$ spherical neighbor coordination shell** (percentile density proxy for relative solvent accessibility) |
| **Hydrophobicity Shift** | $\Delta\text{Hydropathy}$ | **Kyte-Doolittle hydropathy index** (Kyte & Doolittle, *J. Mol. Biol.* 1982) |
| **Steric Volume & Cavities** | $\Delta\text{Volume}$ ($\text{Å}^3$) | **Chothia & Harpaz van der Waals amino acid volumes** (Chothia 1975; Harpaz et al., *Structure* 1994) |
| **Electrostatic Potential** | Net Charge $\Delta$ | **Formal charge at physiological pH 7.0** ($+1$ for Arg/Lys, $-1$ for Asp/Glu) |
| **Thermodynamic Stability** | $\Delta\Delta G$ ($\text{kcal/mol}$) | **PoPMuSiC statistical potentials** (Dehouck et al., *BMC Bioinformatics* 2011) |
| **Sequence Properties** | MW, $\text{pI}$, Instability | **BioPython ProtParam** (Bjellqvist $\text{p}K_a$ values, Guruprasad instability index) |
| **Mechanistic Explanation** | Plain-language summary | **LLM Synthesis (Qwen 2.5 / Gemini Flash)** grounded strictly in computed structural constraints |

---

## Getting Started

### Prerequisites
* Python 3.10 or higher
* A Telegram account (to interact with your bot)

### 1. Clone the Repository
```bash
git clone https://github.com/your-username/fold-bot.git
cd fold-bot
```

### 2. Create and Activate a Virtual Environment
```bash
# Windows
python -m venv .venv
.\.venv\Scripts\activate

# Linux / macOS
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

### 4. Configure Your Environment Variables
Copy the template configuration:
```bash
# Windows PowerShell
Copy-Item env.example .env

# Linux / macOS
cp env.example .env
```

Open `.env` in any text editor and configure your keys:

```env
# 1. Telegram Bot Token (from @BotFather on Telegram)
TELEGRAM_BOT_TOKEN="your_bot_token_here"

# 2. Your Telegram Numeric User ID (message @userinfobot or run /whoami to find it)
ALLOWED_USER_IDS="your_telegram_id"

# 3. LLM Configuration (Choose Gemini or Groq)
LLM_PROVIDER=gemini
GEMINI_API_KEY="your_google_ai_studio_key"
GEMINI_MODEL=gemini-flash-latest

# 4. Automatic Groq Fallback (Free open-source model when Gemini is at capacity)
GROQ_API_KEY="gsk_your_groq_key"
GROQ_MODEL=qwen/qwen3.8-27b

# 5. Structure Backend
FOLD_BACKEND=api
```

>  **Where to get free API keys:**
> * **Telegram Bot Token:** Message [@BotFather](https://t.me/BotFather) on Telegram and type `/newbot`.
> * **Google Gemini API Key:** Free at [aistudio.google.com](https://aistudio.google.com/).
> * **Groq API Key:** Free at [console.groq.com](https://console.groq.com/) (provides lightning-fast Qwen 2.5 inference).
> * **NVIDIA API Key (Optional):** Free at [build.nvidia.com](https://build.nvidia.com/meta/esmfold) for dedicated cloud ESMFold.

### 5. Run the Bot
```bash
python bot.py
```

---

## How to Use FoldBot

You can interact with FoldBot using **slash commands** or **natural plain English**.

### Slash Commands Reference

| Command | Usage | Description |
| :--- | :--- | :--- |
| `/afdb <ID>` | `/afdb P69905` | Fetches precomputed AlphaFold2 structure from EMBL-EBI. |
| `/fold <sequence>` | `/fold MQIFVKTL...` | Predicts 3D structure using ESMFold, returns `.pdb` and cartoon render. |
| `/mutate <protein> <mutation>` | `/mutate P69905 V11R` | Runs structural superposition, Cα RMSD, pLDDT plots, PoPMuSiC $\Delta\Delta G$, and LLM summary. |
| `/stability <protein> <mutation>` | `/stability P69905 V11R` | Computes thermodynamic $\Delta\Delta G$ stability estimates (PoPMuSiC, FoldX, ThermoMPNN). |
| `/analyze <protein or seq>` | `/analyze P69905` | Computes MW, isoelectric point ($\text{pI}$), GRAVY index, and flags low-confidence regions. |
| `/pdb <PDB_ID>` | `/pdb 1UBQ` | Downloads experimental PDB structure from RCSB and renders it. |
| `/method` | `/method` | Opens interactive buttons to switch prediction backends (ESMFold API, NVIDIA, Local). |
| `/llmtest` | `/llmtest` | Tests primary LLM and Groq fallback connections. |
| `/status` | `/status` | Displays active prediction method, compute device, and installed tools. |
| `/reset` | `/reset` | Clears conversation memory and resets the active protein context. |
| `/whoami` | `/whoami` | Shows your Telegram user ID. |

### Natural Language Chatting
Because of the conversational LLM integration, you can also talk to FoldBot naturally:
* *"Fold this sequence: MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"*
* *"What happens if I change residue 11 to Arginine on P69905?"*
* *"How stable is this mutation?"*
* *"What is the molecular weight and isoelectric point of this protein?"*
* *"Explain what pLDDT means in simple terms."*

---

## Security & Privacy

* **Strict Access Control:** Only Telegram user IDs listed in `ALLOWED_USER_IDS` can interact with the bot. Unauthorized users are politely turned away, preventing unauthorized API usage.
* **Never Commit `.env`:** The `.gitignore` is pre-configured to ensure your bot token and API keys are never accidentally pushed to GitHub.
* **Secure Environment Loading:** Environment variables are strictly parsed at startup.

---

## Project Structure

```
fold-bot/
├── bot.py             # Main Telegram bot logic, prediction queues, rendering & analysis pipelines
├── env.example        # Clean environment template with full configuration documentation
├── requirements.txt   # Python package dependencies
├── .gitignore         # Prevents secrets, virtual environments, and cache from being committed
├── LICENSE            # MIT License
└── README.md          # Project documentation and scientific reference
```

---

## Disclaimer & Scientific Limitations

FoldBot is intended for **research, educational, and computational screening purposes only**. 
* Machine learning structure predictors and heuristic statistical potentials generate **hypotheses**, not experimental ground truths.
* Structure predictors are known to have reduced accuracy when ranking minor conformational shifts caused by single point mutations.
* For critical therapeutic or engineering applications, predictions should always be experimentally validated via circular dichroism, differential scanning fluorimetry (DSF), or X-ray crystallography / cryo-EM.

---

## License

This project is licensed under the [MIT License](LICENSE).

