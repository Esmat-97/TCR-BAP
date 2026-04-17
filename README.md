# 🧬 TCR-BAP

**Baseline-Aware Analysis of Antigen-Specific TCR Repertoires**

---

## 📌 Overview

Understanding T-cell receptor (TCR) specificity requires more than identifying frequent patterns — it demands distinguishing true biological signals from baseline repertoire biases.

**TCR-BAP** provides a framework for analyzing antigen-specific TCR repertoires (e.g., CMV) by comparing them against baseline (healthy) repertoires to uncover **true enrichment signals**.

---

## 🎯 Key Idea

Most observed patterns in TCR data (e.g., V gene usage, motifs) may simply reflect underlying repertoire biases.

👉 This project focuses on:

* **Baseline correction**
* **Enrichment analysis**
* **Interpretable specificity profiling**

---

## 🧩 Features

### 1️⃣ TRBV Gene Enrichment

* Computes normalized V gene frequencies
* Calculates **log2 enrichment vs baseline**
* Identifies truly overrepresented genes

### 2️⃣ CDR3 Length Profiling

* Compares length distributions between conditions
* Detects structural constraints in antigen recognition

### 3️⃣ Motif Analysis (Baseline-aware)

* Evaluates motif frequencies relative to baseline
* Avoids misleading interpretations of common motifs

### 4️⃣ Visualization

* Enrichment barplots
* Comparative distributions
* Clean, interpretable figures

---

## 📊 Example Insights

* Selective enrichment of specific TRBV genes (e.g., TRBV29-1)
* Preference for shorter CDR3 lengths (~11–14 aa)
* Motif patterns largely reflect baseline repertoire bias

---

## 🛠️ Tech Stack

* Python 🐍
* Pandas
* NumPy
* Matplotlib / Seaborn



---

## 🚀 How It Works

### Step 1: Load Data

* Antigen-specific TCR dataset (e.g., CMV)
* Baseline repertoire dataset

### Step 2: Normalize Frequencies

Convert raw counts into frequencies

### Step 3: Compute Enrichment

```python
log2_enrichment = log2(CMV_freq / baseline_freq)
```

### Step 4: Visualize Results

* Identify enriched/depleted features
* Interpret biological relevance

---

## 🧠 Why This Matters

Traditional TCR analyses often:
❌ Focus on raw frequencies
❌ Overinterpret common motifs

This framework:
✅ Highlights **true specificity signals**
✅ Provides **interpretable outputs**
✅ Bridges computational analysis with immunological insight

---

## 🔬 Future Directions

* Alpha vs Beta chain comparison
* Cross-reactivity prediction between epitopes
* Integration with structural modeling (e.g., AlphaFold)
* Development of a scoring model (TEMPO-like)

---

## 📎 Data Sources

* VDJdb
* Public TCR repertoire datasets

---

## 🤝 Contributing

Contributions, suggestions, and collaborations are welcome.

---

## 📬 Contact

Mohamed Esmat
Bioinformatics & Immunology

---

## ⭐ Acknowledgment

Inspired by recent advances in interpretable TCR specificity modeling and baseline-aware repertoire analysis.

---
