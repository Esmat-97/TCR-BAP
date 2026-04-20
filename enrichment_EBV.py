import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_EBV = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل gene
human_freq = df_human["TRBV"].value_counts(normalize=True)
EBV_freq = df_EBV["TRBV"].value_counts(normalize=True)

# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": EBV_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["EBV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالـ genes اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("TRBV Gene")
plt.title("TRBV Gene Enrichment in EBV relative to Human baseline")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()






# قراءة البيانات
df_human = pd.read_csv("VDJdb_Human_TRA_TRB_Baseline_Input.csv", sep="\t")
df_EBV = pd.read_csv("VDJdb_EBV_TRA_TRB_Input.csv", sep="\t")

# حساب frequency لكل طول
human_freq = df_human["CDR3b"].str.len().value_counts(normalize=True)
EBV_freq = df_EBV["CDR3b"].str.len().value_counts(normalize=True)

# دمج الاتنين
df = pd.DataFrame({"Human_freq": human_freq, "EBV_freq": EBV_freq}).fillna(0)

# حساب log2 enrichment
df["log2_enrichment"] = np.log2((df["EBV_freq"]+1e-9) / (df["Human_freq"]+1e-9))

# فلترة: نحتفظ بالأطوال اللي Human_freq > 0.01
df = df[df["Human_freq"] > 0.01]

# ترتيب
df = df.sort_values("log2_enrichment", ascending=False)

# رسم barplot
plt.figure(figsize=(12,6))
sns.barplot(x=df.index, y=df["log2_enrichment"], palette="coolwarm")
plt.axhline(0, color="black", linestyle="--")
plt.ylabel("log2 enrichment (EBV vs Human)")
plt.xlabel("CDR3 Length (aa)")
plt.title("CDR3 Length Enrichment in EBV relative to Human baseline")
plt.tight_layout()
plt.show()
