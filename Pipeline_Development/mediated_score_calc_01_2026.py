import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import TunedThresholdClassifierCV, train_test_split
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import average_precision_score
from sklearn import metrics
from sklearn.metrics import r2_score
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import mean_squared_error
import subprocess
import argparse
import os
import math 
from sklearn.metrics import recall_score
from sklearn.metrics import brier_score_loss



parser = argparse.ArgumentParser()
#name of files
parser.add_argument("--samplesheet",help="data_nickname")
parser.add_argument("--snpsheet")
parser.add_argument("--plink_prefix",help="plink files")
parser.add_argument("--skip_raw",action="store_true",help="skip making raw files?")

# whether or not to caulcate, auc and make plots or just return scores
parser.add_argument("--calculate_auc",action="store_true",help="calculate auc & make plots?")
parser.add_argument("--type_of_pheno",help="binary or quantitative phenotype")
parser.add_argument("--use_pcs")
parser.add_argument("--use-covars",help="additional covariates used to test prs",default="")

args = parser.parse_args()

samplesheet_path = args.samplesheet
snpsheet_path = args.snpsheet
plink_prefix = args.plink_prefix
calc_auc = args.calculate_auc
type_of_pheno = args.type_of_pheno
use_pcs = args.use_pcs
skip_raw = args.skip_raw
if args.use_covars != "":
    use_covars = args.use_covars.split(",")
else:
    use_covars = []
print(f'Value of skip_raw: {skip_raw}')



def make_raw_files(plink_input):
    print(plink_input)
    os.environ["plink_prefix"] = plink_input
    os.environ["sample_names"] = "sample_names.txt"
    os.environ["snplist"] = "snp_list.txt"
    os.environ["effect_allele_file"] = f"plink_{samplesheet_path.split('/')[-1]}_ea_file.txt"
    os.environ["plink_output_prefix"] = plink_input.split("/")[-1]
    os.system("bash genotype_filter_01_2026.sh")

def make_log_reg_classifier_overall(input_df, score, pheno, title, covars, pcs):
    if pcs > 0:
        n = pcs + 1
        df_pc_list = list(range(1, n))
        df_pc_formatted_list = [f"PC{n}" for n in df_pc_list]
        df = input_df[[score] + [pheno] + df_pc_formatted_list + covars].copy()
        df[score] = (df[score] - df[score].mean()) / df[score].std()
        df = df.dropna()
        df = df.astype('float64')
        print(df.columns)
        X = df[[score] + df_pc_formatted_list + covars]
        y = df[pheno]
    else:
        df = input_df[[score] + [pheno]].copy()
        df[score] = (df[score] - df[score].mean()) / df[score].std()
        df = df.dropna()
        df = df.astype('float64')
        print(df.columns)
        X = df[[score]]
        y = df[pheno]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25,random_state=0)
    clf = LogisticRegression(random_state=0,multi_class='auto',solver='liblinear')
    clf.fit(X_train, y_train)
    clf_tuned = TunedThresholdClassifierCV(
        clf, scoring="balanced_accuracy"
        ).fit(X_train, y_train)
    y_preds = clf_tuned.predict_proba(X_test)
    y_probs = y_preds[:,1]
    y_pred_prob_binary = (y_probs > 0.5).astype(int)

    #AUROC
    fpr, tpr, thresholds = metrics.roc_curve(y_test,y_probs)
    roc_auc = metrics.auc(fpr, tpr)
    print(roc_auc)
    y_test_vals = y_test.values.ravel()
    n1 = sum(y_test_vals == 1)
    n2 = sum(y_test_vals != 1)
    q1 = roc_auc / (2 - roc_auc)
    q2 = 2 * roc_auc**2 / (1 + roc_auc)
    se_roc_auc = math.sqrt((roc_auc*(1 - roc_auc) + (n1 - 1)*(q1 - roc_auc**2) + (n2 - 1)*(q2 - roc_auc**2)) / (n1*n2))
    roc_auc_ci = 1.96*se_roc_auc
    #PRAUC
    prauc =  average_precision_score(y_test,y_probs)
    # Accuracy
    accuracy = clf.score(X_test,y_test)
    n = len(y_test)  # Total number of samples
    correct = int(accuracy * n)  # Number of correct predictions
    alpha = 0.05  # for 95% confidence interval
    z = 1.96  # z-score for 95% confidence
    denominator = 1 + z**2/n
    p_adj = (correct + z**2/2)/denominator
    spread = z * math.sqrt(max(0, (accuracy * (1 - accuracy) + z**2 / (4 * n)) / denominator))
    acc_conf_int = spread

    #Balanced Accuracy

    balanced_accuracy = balanced_accuracy_score(y_test,y_pred_prob_binary)
    correct_balanced = int(balanced_accuracy * n)  # Number of correct predictions
    denominator_balanced = 1 + z**2/n
    p_adj_balanced = (correct_balanced + z**2/2)/denominator_balanced
    balanced_accuracy_spread = z * math.sqrt(max(0, (balanced_accuracy * (1 - balanced_accuracy) + z**2 / (4 * n)) / denominator_balanced))
    bal_acc_conf_int = balanced_accuracy_spread


    #sensity & specificity calc
    sensitivity = recall_score(y_test, y_pred_prob_binary)
    specificity = recall_score(np.logical_not(y_test), np.logical_not(y_pred_prob_binary))
    #brier score
    brier_score = brier_score_loss(y_test, y_pred_prob_binary)

    print(f"========={title} Performance Metrics=========")
    print(f"{title} AUROC:{roc_auc:.4f }\u00b1{roc_auc_ci:.4f}")
    print(f"{title} Accuracy:{accuracy:.4f}\u00b1{acc_conf_int:.4f}")
    print(f"{title} Balanced Accuracy: {balanced_accuracy:.4f}\u00b1{bal_acc_conf_int:.4f}")
    print(f"{title} PRAUC:{prauc:.4f}")
    print(f"{title} Sensitivity:{sensitivity:.4f}")
    print(f"{title} Specificity:{specificity:.4f}")
    print(f"{title} brier_score:{brier_score:.4f}")

    performance_metrics = {"AUROC":[roc_auc],"Accuracy":[accuracy],"Sensitivity":[sensitivity],
    "Specificity":[specificity],"brier_score":[brier_score],"PRAUC":[prauc]}
    performance_metrics_df = pd.DataFrame.from_dict(performance_metrics)

    display = RocCurveDisplay.from_estimator(clf_tuned, X_test, y_test,name=f'{score},title:{title}')
    display.plot()
    plt.savefig(f"./{title}_AUC_ROC.png")  
    return clf,roc_auc,display,performance_metrics_df


def make_lin_reg_classifier_overall(input_df,score,pheno,title,covars,pcs):
    if pcs > 0:
        n = pcs + 1
        df_pc_list = list(range(1, n))
        df_pc_formatted_list = [f"PC{n}" for n in df_pc_list]
        df = input_df[[score] + [pheno] + df_pc_formatted_list + covars].copy()
        df[score] = (df[score] - df[score].mean()) / df[score].std()
        df = df.dropna()
        df = df.astype('float64')
        print(df.columns)
        X = df[[score] + df_pc_formatted_list + covars]
        y = df[pheno]
    else:
        df = input_df[[score] + [pheno]].copy()
        df[score] = (df[score] - df[score].mean()) / df[score].std()
        df = df.dropna()
        df = df.astype('float64')
        print(df.columns)
        X = df[[score]]
        y = df[pheno]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25,random_state=0)
    clf = LinearRegression(positive=True)
    y_preds = clf.fit(X_train,y_train).predict(X_test)

    #R^2
    rscore = r2_score(y_test,y_preds)
    se = sqrt((4*rscore*((1-rscore)**2)*((n-k-1)**2))/((n**2-1)*(n + 3)))
    r2ci = 2*se

    #MSE
    mse = mean_squared_error(y_test, y_preds)

    print(f'========={title} Performance Metrics=========')
    print(f"R^2: {rscore:.4f}\u00b1{r2ci:.4f}")
    print(f"Mean Squared Error:{mse:.4f}")

    plt.plot(X_test, y_preds, color="blue", linewidth=3)
    plt.xticks(())
    plt.yticks(())
    display = plt.show()
    return clf,rscore,display




# make sample names file for plink input
samplesheet  = pd.read_table(samplesheet_path,index_col="IID")

fam_file = pd.read_table(f"{plink_prefix}.fam",usecols=[0,1],header=None,names=["FID","IID"],sep="\s+")
sample_names = fam_file
sample_names.to_csv("./sample_names.txt",sep="\t",index=False)


#make snp_list file for plink_input
print("Reading SNP List...")
snp_sheet = pd.read_table(snpsheet_path,sep="\s+")
print(snp_sheet.columns)
snp_list = pd.DataFrame()
snp_list["SNP_ID"] = snp_sheet["SNP_ID"]
print(snp_list.columns)
print(snp_list)
if snp_list["SNP_ID"].str.contains(":").any():
    snp_list[["chr", "start"]] = snp_list["SNP_ID"].str.split(":").str[:2].tolist()
elif snp_list["SNP_ID"].str.contains("_").any():
    snp_list[["chr", "start"]] = snp_list["SNP_ID"].str.split("_").str[:2].tolist()
snp_list["stop"] = snp_list["start"]
print("LOOK HERE!")
print(snp_list.columns)
#snp_list['SNP_ID'] = snp_list['chr'].str.replace('chr', '', regex=False) + ':' + snp_list['start'].astype(str)
snp_ids = snp_list["SNP_ID"]
snp_toextract = snp_list['SNP_ID']
snp_toextract.to_csv("snp_list.txt",sep="\t",index=False,header=False)
effect_allele_file = pd.DataFrame()
effect_allele_file['SNP_ID'] = snp_sheet['SNP_ID']
effect_allele_file['A1'] = snp_sheet['EA_MP']
effect_allele_file.to_csv(f"plink_{samplesheet_path.split('/')[-1]}_ea_file.txt",sep="\t",index=False)

# make raw files to use as genotype df
if not skip_raw:
    print("Making Raw Files...")
    make_raw_files(plink_prefix)
else:
    print("Using Existing Raw Files")

plink_output_prefix = plink_prefix.split("/")[-1]

filtered_genotypes = pd.read_csv(f"{plink_output_prefix}_filtered.raw",index_col='IID',sep="\s+" )

print(f'this is the filtered_geno_df {filtered_genotypes.head()}')

#print(filtered_genotypes.columns)
# process 1 parse input sheet and split into two groups to create dictionaries of diff snps
mediated_phenotype_snps = {}
mediator_snps = {}
print("SNP IDS CHECK")
print(snp_sheet["SNP_ID"].head(10))
# dictionaries with SNPIDS formatted to match plink 2 output format
mediated_phenotype_snps.update({
    row["SNP_ID"] + (f"_{row['EA_MP']}"): row["BETA_MP"]
    for _, row in snp_sheet.iterrows()
}) # maintain betas from mediated phenotyype GWAS

mediator_snps.update({
    row["SNP_ID"] + (f"_{row['EA_M']}"): row["BETA_M"] if row["EA_M"] == row["EA_MP"] else -row["BETA_M"]
    for _, row in snp_sheet.iterrows() # if the EA_M is the same as EA_MP then the Mediator GWAS Beta is maintained if not it is reverse
})


duplicates = filtered_genotypes.index[filtered_genotypes.index.duplicated()]
# Print duplicates
print("Filtered Genotypes Duplicate indices:", duplicates.tolist())
duplicates = samplesheet.index[samplesheet.index.duplicated()]
# Print duplicates
print("Samplesheeet Duplicate indices:", duplicates.tolist())
samplesheet = samplesheet[~samplesheet.index.isna()]
print(samplesheet.index.isna())


#make df with genotype and mediator 
genotype_and_mediator = pd.concat([filtered_genotypes,samplesheet],axis=1,join="inner")
print(f"This is the dataframe used to calculate PRS")
print(f"{genotype_and_mediator.head()}")

formatted_snp_ids = snp_sheet["SNP_ID"] + "_" + snp_sheet["EA_MP"] # add effect allele to snp ids to match with raw file snp ids
print("THESE ARE FORMATTED SNP IDS")
print(formatted_snp_ids.head(10))


print("Making File for PRS Calculation")
mediated_phenotype_col_list = []
mediator_col_list  = []

clean_snp_list = [snp for snp in formatted_snp_ids if snp in list(genotype_and_mediator.columns)]
dropped_snps = list(set(formatted_snp_ids) - set(clean_snp_list))
dropped_snps = pd.Series(data=dropped_snps)


print(f"{(len(clean_snp_list)/ len(formatted_snp_ids)) * 100} % of variants matched")
print(f"{(dropped_snps.shape[0]/len(formatted_snp_ids))* 100} % of variants dropped")

for snp in clean_snp_list:
    genotype_and_mediator[f"{snp}_m"] = genotype_and_mediator[f"{snp}"] * mediator_snps[f"{snp}"] # multiply genotype by mediator beta
    genotype_and_mediator[f"{snp}_mp"] = genotype_and_mediator[f"{snp}"] * mediated_phenotype_snps[f"{snp}"] # multiply genotype by mediated phenotype beta
    if mediator_snps[f"{snp}"]!=0:
        mediator_col_list.append(snp)
    if mediated_phenotype_snps[f"{snp}"] != 0:
        mediated_phenotype_col_list.append(snp)

dropped_snps.to_csv("Dropped_SNPs.txt",index=False)
prs_values = genotype_and_mediator[["M","P"]].copy()
print(genotype_and_mediator.columns)

print("STARTING PRS VALUE DF")

use_pcs = int(use_pcs)

# PRS_P = BETA_MP PRS
prs_values['PRS_P'] = genotype_and_mediator[mediated_phenotype_col_list].sum(axis=1) 
# PRS_M = Mediated Phenotype Individual Value * BETA_M PRS
prs_values['PRS_M'] = genotype_and_mediator["M"].values * genotype_and_mediator[mediator_col_list].sum(axis=1)
# Mediated PRS sum of PRS_M for PRS_P
prs_values["PRS"] = prs_values[['PRS_P', 'PRS_M']].sum(axis=1)
if use_pcs > 0:
    prs_values[[f"PC{n}" for n in range(1, use_pcs + 1)]] = genotype_and_mediator[[f"PC{n}" for n in range(1, use_pcs + 1) ]].copy()
if len(use_covars) > 0: 
    prs_values[use_covars] = genotype_and_mediator[use_covars].copy()
prs_values = prs_values.dropna()

print(f"PRS Values Caluclated for {prs_values.shape[0]} individuals")

folder = "./"
filename = f"{snpsheet_path.split('/')[-1].replace('txt','')}_Mediated_PRS_output.txt"
full_path = os.path.join(folder, filename)

prs_values.to_csv(full_path, sep="\t")

print(f"PRS Values Caluclated saved to {full_path}")

if calc_auc and not prs_values.empty:
    if type_of_pheno == "binary":
        print("PRS Calculated Generating Plots & Calculating AUC")
        samplesheet_title = samplesheet_path.split("/")[-1].replace(".txt","")
        # output = make_log_reg_classifier_overall(prs_values,"PRS","P",samplesheet_title)
        clf_prs,roc_auc_prs,display_prs,performance_metrics_prs = make_log_reg_classifier_overall(prs_values,"PRS","P",samplesheet_title,use_covars,use_pcs)
        performance_metrics_prs.to_csv(f"{samplesheet_title}_metrics.txt",index=False,sep="\t")
    elif type_of_pheno == "quantitative":
        print("PRS Calculated Generating Plots & Calculating R^2")
        samplesheet_title = samplesheet_path.split("/")[-1].replace(".txt","")
        output = make_lin_reg_classifier_overall(prs_values,"PRS","P",samplesheet_title,use_covars,use_pcs)
elif prs_values.empty:
    print("PRS Table Empty Check Execution")   
elif not calc_auc and not prs_values.empty:
    print(f"PRS Calculated and printed to {samplesheet_path}_{snpsheet_path}_Mediated_PRS.txt")
else:
    print("Bad Selection.Check execution")