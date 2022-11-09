## Introduction

-   

[1]

-   

[2-7]

Upon stimulation, B cells switch from naïve B cells with IgM or IgD antibodies/receptors to (eventually) antibody-producing plasma cells secreting IgA, IgG, or IgE. Whereas *affinity* for an antigen is improved through somatic hypermutation in the variable region (mostly in the V and J genes), this process of changing Ig type, termed class switching, is achieved through recombination in the heavychain constant region, and improves a receptor’s *effector* functions.[8, 9]

## 

## 

## 

## Materials and Methods

### Sample collection

A total of peripheral blood mononuclear cells (PBMCs) were analysed: samples were collected from 6 COVID-19 inpatients admitted into Juntendo University Hospital in Tokyo, Japan, at several points in time with different stages of disease progression (see **Table S1** for patient information); and 3 samples were collected from 3 healthy donors. The median age for patients was 76 years (range: 53-93), and all but Patient 1 were male.

As shown in Table S1, P1 to 4 are cases in which specimens were collected at the time of transition from moderate to critical stage, and P5 and 6 are cases in which specimens were collected at the stage of transition from critical to moderate. P1, 2, and 3 are dead cases, and Pt4, 5, and 6 are surviving cases.

SARS-CoV-2 infection was established through RT-PCR-based molecular testing,and was performed using nasopharyngeal specimens by the 2019 Novel Coronavirus Detection Kit (Shimadzu, Kyoto, Japan).[10]

[11][11]

### Sample processing

PBMCs were isolated using BD Vacutainer CPT tubes according to the manufacturer’s recommendations. The collected cells were suspended in 1 mL of CELLBANKER 1 plus (Nihon Zenyaku Kogyo, Koriyama, Japan), frozen at -80°C, and then used for flow cytometry analysis and 10x single-cell RNA sequencing (scRNA-seq). Serum samples were stored at -80°C prior to use.

### Single cell RNA-seq assay (10x platform)

PBMC suspensions were barcoded through the 10x Chromium Single Cell platform using Chromium Next GEM single cell 5’ GEM reagent kits v2 (PN-1000263, 10x Genomics, CA) according to the manufacturer’s protocol. The loaded cells numbered approximately 1.6 × 104, aiming for 1.0 × 104 single cells per reaction on the Chip K Single Cell Kit (PN-1000286, 10x Genomics) with the Chromium controller (10x Genomics). Single-cell RNA libraries were prepared, and unique sample indexes added using the library construction kit, the Single Index Kit T Set A (PN-1000213, 10x Genomics), and a universal library conversion kit (App-A, cat. No. 1000004155, MGI, China). The libraries were sequenced using the DNBSEQ-G400 platform (MGI).

### Data processing

Raw sequencing reads in FASTQ format were analyzed using the Cell Ranger (10x Genomics, version 6.1.2) multi pipeline, performing gene expression (GEX) alignment against the GRCh38 human reference dataset (10x Genomics) with a SARS-CoV-2 genomic RNA (GenBank: LC606022.1) appended. V(D)J assembly and alignment and paired clonotype calling for T and B cell receptors were done against the refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 reference (10x Genomics). CDR3 sequences and rearranged full-length BCR/TCR V(D)J segments as well as clonotype frequency were also obtained.

### Single-cell RNA-seq

Downstream single-cell data analysis of Cell Ranger’s output expression count matrices was performed using the R programming language (version 4.2.1)[12] with the Seurat package (version 4.0).[13] All samples were integrated into one dataset, categorized according to severity and patient identity. Cell cycle phases were determined for each cell, which were to have no discernible effect on the variation of the data. Upon evaluation of the quality of the gene counts, the collective sample was restricted to cells associated with more than 500 transcripts, more than 250 genes, and less than 20% mitochondrial genes. Genes with zero counts were also removed from the dataset. After this filtration, 19,9 cells were kept for further analysis.

The data were normalized, and variable features determined for clustering the cells. The Azimuth web application was used for cluster annotation,[13] and UMAP dimensionality reduction was used to view the resulting cell clusters/types. Azimuth annotations were used as cell type identifiers in further analysis.

### Immune repertoire analysis

Gene expression data matrices from the scRNA-seq analysis were subset to two groups: T cells (CD3-enriched cells), and B cells. To these were integrated the clonotype filtered contig data for each library in the scRepertoire R package (version 1.7.0)[14] using the *CombineTCR* and *CombineBCR* functions. Only productive TCR/BCR sequences were kept for subsequent analysis. TCR/BCR repertoire analysis was performed using the Seurat and scRepertoire packages in R.

### 

[15][16][17][18][19][20]

### 

[21][22]

[23][24][25][11]

### Statistical analysis

The Shannon index was used to calculate the V gene usage and Ig isotype diversity scores for the individual samples[26] and Wilcoxon’s test was used to compare the clonal and isotype diversity indices between samples, to compare isotype proportions, and to detect differentially expressed genes (DEGs) between different subsets of the data.[27] In all statistical comparisons, *p* values less than 0.05 were considered significant. All statistical analyses were performed using the R programming language.[12]

## Results

Thirty cell clusters were identified by the Azimuth application across 19, cells. Of these, 30,071 (2.%) were from the samples of the healthy donors, 15,493 (%) were from the COVID-19 patients’ samples in mild stage, 32,182 (2%) in moderate stage, and 51,723 (%) in critical stage (**Figure 1**). Fine resolution clustering of mRNA profiles revealed 13 clusters of T lymphocytes that include CD4 and Treg subtypes, CD8 subtypes, MAIT cells, double-negative T cells (dnT), and gamma-delta T cells (gdT)

### 

### 

### 

[28]

[29]

### 

V-J usage pattern for CD4 TCM cells expanded common clonotypes for infected patients and even in controls, as well as a more diverse usage combination of the V-gene with the J region. But CD8 TEM V-J usage was extremely variable for each patient, and V-J usage was less diverse.

Taken together,

### 

### 

### 

ethe has been analyzed in enriched T cells We focused on CD4 TCM cells which ut.

To investigate While all patients underwent critical state around 11-14 days from onset of the COVID-19 infection, the moderate/mild stages were of two different timepoints. Patient 1–4 were sampled during the early stages of the disease (2-5 days from onset) before progressing to the critical state. Patient 5 and 6, however were sampled during the recovery phase of the disease (41 days and 29 days, respectively). In this section, we performed differential gene expression analysis with comparison between the progressing group vs the recovering group to determine upregulated genes and pathways for each group.

along with in moderate stage compared to critical stageWe further performed the to understand the ( observed that findings were consistent for all s. in almost all patients

### 

### 

identifiedThen

two were detectedIn thecases, andwere activated

### B

Integration of the V(D)J sequences into the scRNA-seq data 6 cellsPlasmablasts (33.%), B naïve (35.%), B memory (9.6%), and B intermediate (21.%).

t all stages of COVID-19 patients examined in this study, the percentage of plasmablasts was significantly higher than in healthy controls. However, changes in the proportion of plasmablasts according to severity differed from case to case and no consistent trends were observed. Pt14 progressed from moderate to critical stages, Pt1 2 increased as they progressed to critical stages, in Pt3 and, the plasmablast ratio at critical stages decreased compared to moderate stages. In Pt5 and 6, who recovered from critical to moderate stage, plasmablasts increased and decreased during their recovery process, respectively. Of note, intermediate B cells showed changes opposite to increases and decreases in plasmablasts in all cases. The proportion of B naive cells decreased with progression from moderate to critical in dead cases (Pt1, Pt2, Pt3), but increased in survivor (Pt4). In survivors Pt5 and Pt6, the proportion of B naive cells decreased with recovery from critical to moderate stage.

### 

BCR V gene usage diversity (p value = 0.0029). The proportion of expanded clones (3 cells or more) was significantly higher in critical cases than in moderate ones (p =0.008).

The proportion of expanded clones (3 cells or more) was significantly higher in critical/severe cases than in mild/moderate ones (p =0.008).

The most abundant V–J gene usage combination was IGHV3-23–IGHJ4, which was in all samples regardless of disease status and severity. The second highest frequency was with IGHV4-34–IGHJ5, mostly detected in the deceased Patient 2’s critical sample (**Figure 2**). Of particular interest was the IGH1-18 gene, most commonly with IGHJ3 in the deceased Patient 1 in the critical (and terminal) stage. Also of interest was IGHV4-39, most highly used in Patient 3’s (deceased) critical sample. This gene has also been used for CLL prognostication, and as a marker for progression into Richter Syndrome.[30]

J gene usage analysis revealed preferential usage of the IGHJ4 gene for all samples, except for Patient 2, whose critical sample showed a predominance of IGHJ5 usage. No significant difference was in J gene usage diversity across the samples.

Concerning isotype diversity, more diverse isotypes in COVID-19 patients than in healthy controls (p = 0.0059), which is the opposite direction to V gene usage diversity. No significant difference in isotype diversity was found between deceased and recovered patients. In our collective data set, the most common isotype was IgM; it was also the most common Ig type in the individual samples except those of Patient 2 and Patient 3—both deceased. Patient 2 had larger proportions of IgA1, IgA2, and IgG1; and Patient 3 had IgA1 in the largest proportion. The use of IgA1 was significantly higher in the deceased than in the recovered patients (p = 0.0013) and healthy controls (p = 0.024)

### Differential Gene Expression

### 

[11, 12][13][14]

## 

## 

## 

## 

It has been suggested that Covid-19 patients exhibit distinct immune responses to SARS-CoV-2 infection with varying severities.[31] In this study, we sought to better understand differences in the immune response among different severity groups, as well as patients with different outcomes.

-   

[32][3][7]

-   

[33, 34]

[35][3, 35, 36][37]

-   

[38][39]

-   

T cells in general have been described to sustain inflammation with elevated interferon (INF)-α, tumor necrosis factor (TNF), and INF-γwith increasing SARS-CoV-2 severity.

[40-42]

-   

It has been suggested that Covid-19 patients exhibit distinct immune responses to SARS-CoV-2 infection with varying severities.[34] In this study, we sought to better understand differences in the immune response among different severity groups, as well as patients with different outcomes.

[44][31][45]

In B cell repertoire analysis, special emphasis was given to the most common heavy-chain Immunoglobulin Heavy Variable (IGHV) genes, especially in combination with J genes; as these are thought to drive the most variation in B cell receptors.[46] The most common V gene in our data set, IGHV3-23, has been reported in other COVID-19 studies as well.[44, 47, 48] IGHV3-23 has also been found to be associated with aggressive chronic lymphocytic leukaemia (CLL) where it is highly mutated,[49] indicating bad prognosis.[50]

IGHV4-34, the second most commonly used V gene, mostly in one of our deceased patients, is associated with self-reactive antibodies,[51] and has also been previously observed in COVID-19 patients.[31, 48, 52] Interestingly, Galson et al.[52] reported the usage of this V gene to be significantly higher in B cells with the IgG1 isotype in improving COVID-19 patients than those in stable or deteriorating patients, while in our data set, IGHV-34 associated with IgG1 was mostly detected in one of the deteriorating (and eventually deceased) patients.

IGHV3-33 was the third most common V gene, in combination with IGHJ4. Other studies also identified an increase in this gene in association with COVID-19.[44, 52, 53] It is also noteworthy that IGHV3-33 has been detected in potent antibodies in the recipients of malaria vaccination. [54]

IGHV1-18–IGHJ3 was the most common gene segment combination in Patient 1’s critical sample. It has been reported that use of this gene is elevated in acute dengue “without warning signs”[55], and has been identified as first in a panel of biomarkers for predicting complete remission to treatment in patients with proliferative lupus nephritis.[56] IGHV1-18 is also associated with tumor immunity, and others have detected it as a positive prognostic marker gene in immune cells in the melanoma tumor microenvironment.[57, 58]

[59][60][60][61, 62][63, 64]he MYC Targets V1 pathway. This pathway is thought to promote cell proliferation through involvement in mitochondrial biogenesis and oxidative phosphorylation.[65, 66] Enrichment of this pathway in the recovering patients might be indicative of tissue repair after inflammation.

Oxidative hosphorylation. This process is important in innate immunity and T cell proliferation.[67, 68] Dysfunction in this pathway is thought to be associated with excessive release of reactive oxygen species and other inflammatory cytokines and enzymes involved in promoting hyperinflammation in sepsis, and might increase viral replication in diseased cells.[69] It is also believed that, through manipulating the metabolic and signalling functions of mitochondria, SARS-CoVproteins can evade the immune system and alter the host’s anti-viral response.[70] It is possible that this kind of viral “hijacking” of the mitochondria is responsible for the high enrichment of this pathway.

[71-74][75][76, 77]

[78]

### 

### 

## References

1\. Consortium CMBA. A blood atlas of COVID-19 defines hallmarks of disease severity and specificity. Cell. 2022;185(5):916-38.e58. Epub 20220121. doi: 10.1016/j.cell.2022.01.012. PubMed PMID: 35216673; PubMed Central PMCID: PMCPMC8776501.

2\. Chen Z, John Wherry E. T cell responses in patients with COVID-19. Nature Reviews Immunology. 2020;20(9):529-36. doi: 10.1038/s41577-020-0402-6.

3\. Grifoni A, Weiskopf D, Ramirez SI, Mateus J, Dan JM, Moderbacher CR, et al. Targets of T Cell Responses to SARS-CoV-2 Coronavirus in Humans with COVID-19 Disease and Unexposed Individuals. Cell. 2020;181(7):1489-501.e15. Epub 20200520. doi: 10.1016/j.cell.2020.05.015. PubMed PMID: 32473127; PubMed Central PMCID: PMCPMC7237901.

4\. Sekine T, Perez-Potti A, Rivera-Ballesteros O, Strålin K, Gorin JB, Olsson A, et al. Robust T Cell Immunity in Convalescent Individuals with Asymptomatic or Mild COVID-19. Cell. 2020;183(1):158-68.e14. Epub 20200814. doi: 10.1016/j.cell.2020.08.017. PubMed PMID: 32979941; PubMed Central PMCID: PMCPMC7427556.

5\. Mazzoni A, Vanni A, Spinicci M, Capone M, Lamacchia G, Salvati L, et al. SARS-CoV-2 Spike-Specific CD4+ T Cell Response Is Conserved Against Variants of Concern, Including Omicron. Frontiers in Immunology. 2022;13.

6\. Le Bert N, Tan AT, Kunasegaran K, Tham CYL, Hafezi M, Chia A, et al. SARS-CoV-2-specific T cell immunity in cases of COVID-19 and SARS, and uninfected controls. Nature. 2020;584(7821):457-62. Epub 20200715. doi: 10.1038/s41586-020-2550-z. PubMed PMID: 32668444.

7\. Jafari M, Kolahdooz H, Mahmoudi M, Azarnaminy AF, Mobasheri L, Esmaeili SA. The impact of lymphoid memory cells in different ages of COVID-19 patients. Gene Rep. 2022;26:101503. Epub 20220115. doi: 10.1016/j.genrep.2022.101503. PubMed PMID: 35071822; PubMed Central PMCID: PMCPMC8760633.

8\. Alam R, Gorska M. 3. Lymphocytes. J Allergy Clin Immunol. 2003;111(2 Suppl):S476-85. doi: 10.1067/mai.2003.121. PubMed PMID: 12592294.

9\. Alberts BWJHT. Molecular biology of the cell. New York: Garland Science; 2008.

10\. Fukumoto T, Iwasaki S, Fujisawa S, Hayasaka K, Sato K, Oguri S, et al. Efficacy of a novel SARS-CoV-2 detection kit without RNA extraction and purification. Int J Infect Dis. 2020;98:16-7. Epub 20200626. doi: 10.1016/j.ijid.2020.06.074. PubMed PMID: 32599282; PubMed Central PMCID: PMCPMC7318955.

11\. WHO. Living guidance for clinical management of COVID-19 2021. Available from: <https://www.who.int/publications/i/item/WHO-2019-nCoV-clinical-2021-2>.

12\. R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical

Computing. 4.2.1 ed. Vienna, Austria2022.

13\. Hao Y, Hao S, Andersen-Nissen E, Mauck WM, 3rd, Zheng S, Butler A, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573-87.e29. Epub 20210531. doi: 10.1016/j.cell.2021.04.048. PubMed PMID: 34062119; PubMed Central PMCID: PMCPMC8238499.

14\. Borcherding N. scRepertoire: A toolkit for single-cell immune receptor profiling. 1.7.0 ed2022.

15\. Tickotsky N, Sagiv T, Prilusky J, Shifrut E, Friedman N. McPAS-TCR: a manually curated catalogue of pathology-associated T cell receptor sequences. Bioinformatics. 2017;33(18):2924-9. doi: 10.1093/bioinformatics/btx286. PubMed PMID: 28481982.

16\. Shugay M, Bagaev DV, Zvyagin IV, Vroomans RM, Crawford JC, Dolton G, et al. VDJdb: a curated database of T-cell receptor sequences with known antigen specificity. Nucleic Acids Res. 2018;46(D1):D419-d27. doi: 10.1093/nar/gkx760. PubMed PMID: 28977646; PubMed Central PMCID: PMCPMC5753233.

17\. Gielis S, Moris P, Bittremieux W, De Neuter N, Ogunjimi B, Laukens K, et al. Detection of Enriched T Cell Epitope Specificity in Full T Cell Receptor Sequence Repertoires. Frontiers in Immunology. 2019;10. doi: 10.3389/fimmu.2019.02820.

18\. Chronister WD, Crinklaw A, Mahajan S, Vita R, Koşaloğlu-Yalçın Z, Yan Z, et al. TCRMatch: Predicting T-Cell Receptor Specificity Based on Sequence Similarity to Previously Characterized Receptors. Frontiers in Immunology. 2021;12.

19\. Gantner P, Pagliuzza A, Pardons M, Ramgopal M, Routy JP, Fromentin R, et al. Single-cell TCR sequencing reveals phenotypically diverse clonally expanded cells harboring inducible HIV proviruses during ART. Nat Commun. 2020;11(1):4089. Epub 20200814. doi: 10.1038/s41467-020-17898-8. PubMed PMID: 32796830; PubMed Central PMCID: PMCPMC7427996.

20\. Meysman P, De Neuter N, Gielis S, Bui Thi D, Ogunjimi B, Laukens K. On the viability of unsupervised T-cell receptor sequence clustering for epitope preference. Bioinformatics. 2019;35(9):1461-8. doi: 10.1093/bioinformatics/bty821. PubMed PMID: 30247624.

21\. Browning BL. PRESTO: Rapid calculation of order statistic distributions and multiple-testing adjusted P-values via permutation for one and two-stage genetic association studies. BMC Bioinformatics. 2008;9(1):309. doi: 10.1186/1471-2105-9-309.

22\. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417-25. doi: 10.1016/j.cels.2015.12.004. PubMed PMID: 26771021; PubMed Central PMCID: PMCPMC4707969.

23\. Raudvere U, Kolberg L, Kuzmin I, Arak T, Adler P, Peterson H, et al. g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Research. 2019;47(W1):W191-W8. doi: 10.1093/nar/gkz369.

24\. Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Res. 2003;13(11):2498-504. doi: 10.1101/gr.1239303. PubMed PMID: 14597658; PubMed Central PMCID: PMCPMC403769.

25\. Bader GD, Hogue CW. An automated method for finding molecular complexes in large protein interaction networks. BMC Bioinformatics. 2003;4:2. Epub 20030113. doi: 10.1186/1471-2105-4-2. PubMed PMID: 12525261; PubMed Central PMCID: PMCPMC149346.

26\. Keylock CJ. Simpson diversity and the Shannon–Wiener index as special cases of a generalized entropy. Oikos. 2005;109(1):203-7. doi: <https://doi.org/10.1111/j.0030-1299.2005.13735.x>.

27\. Wilcoxon F. Individual comparisons of grouped data by ranking methods. J Econ Entomol. 1946;39:269. doi: 10.1093/jee/39.2.269. PubMed PMID: 20983181.

28\. Quiros-Fernandez I, Poorebrahim M, Fakhr E, Cid-Arregui A. Immunogenic T cell epitopes of SARS-CoV-2 are recognized by circulating memory and naive CD8 T cells of unexposed individuals. EBioMedicine. 2021;72:103610. Epub 20211006. doi: 10.1016/j.ebiom.2021.103610. PubMed PMID: 34627082; PubMed Central PMCID: PMCPMC8493415.

29\. Shrock E, Fujimura E, Kula T, Timms RT, Lee IH, Leng Y, et al. Viral epitope profiling of COVID-19 patients reveals cross-reactivity and correlates of severity. Science. 2020;370(6520). doi: 10.1126/science.abd4250. PubMed PMID: 32994364.

30\. Rossi D, Spina V, Cerri M, Deambrogi C, De Paoli L, Laurenti L, et al. Usage of IGHV4-39 with Stereotypic B Cell Receptor Is An Independent Risk Factor of Chronic Lymphocytic Leukemia Transformation to Richter Syndrome. Blood. 2008;112(11):778. doi: <https://doi.org/10.1182/blood.V112.11.778.778>.

31\. Zhang F, Gan R, Zhen Z, Hu X, Li X, Zhou F, et al. Correction to: Adaptive immune responses to SARS-CoV-2 infection in severe versus mild individuals. Signal Transduct Target Ther. 2021;6(1):161. Epub 20210419. doi: 10.1038/s41392-021-00540-4. PubMed PMID: 33875639; PubMed Central PMCID: PMCPMC8054129.

32\. Gray JI, Westerhof LM, MacLeod MKL. The roles of resident, central and effector memory CD4 T-cells in protective immunity following infection or vaccination. Immunology. 2018. Epub 20180323. doi: 10.1111/imm.12929. PubMed PMID: 29570776; PubMed Central PMCID: PMCPMC6050220.

33\. Sewell AK. Why must T cells be cross-reactive? Nature Reviews Immunology. 2012;12(9):669-77. doi: 10.1038/nri3279.

34\. Kundu R, Narean JS, Wang L, Fenn J, Pillay T, Fernandez ND, et al. Cross-reactive memory T cells associate with protection against SARS-CoV-2 infection in COVID-19 contacts. Nature Communications. 2022;13(1):80. doi: 10.1038/s41467-021-27674-x.

35\. Mahajan S, Kode V, Bhojak K, Karunakaran C, Lee K, Manoharan M, et al. Immunodominant T-cell epitopes from the SARS-CoV-2 spike antigen reveal robust pre-existing T-cell immunity in unexposed individuals. Scientific Reports. 2021;11(1):13164. doi: 10.1038/s41598-021-92521-4.

36\. Simnica D, Schultheiß C, Mohme M, Paschold L, Willscher E, Fitzek A, et al. Landscape of T-cell repertoires with public COVID-19-associated T-cell receptors in pre-pandemic risk cohorts. Clin Transl Immunology. 2021;10(9):e1340. Epub 20210828. doi: 10.1002/cti2.1340. PubMed PMID: 34484739; PubMed Central PMCID: PMCPMC8401425.

37\. Ng OW, Chia A, Tan AT, Jadi RS, Leong HN, Bertoletti A, et al. Memory T cell responses targeting the SARS coronavirus persist up to 11 years post-infection. Vaccine. 2016;34(17):2008-14. Epub 20160305. doi: 10.1016/j.vaccine.2016.02.063. PubMed PMID: 26954467; PubMed Central PMCID: PMCPMC7115611.

38\. Wang P, Jin X, Zhou W, Luo M, Xu Z, Xu C, et al. Comprehensive analysis of TCR repertoire in COVID-19 using single cell sequencing. Genomics. 2021;113(2):456-62. Epub 20201228. doi: 10.1016/j.ygeno.2020.12.036. PubMed PMID: 33383142; PubMed Central PMCID: PMCPMC7833309.

39\. Mallajosyula V, Ganjavi C, Chakraborty S, McSween AM, Pavlovitch-Bedzyk AJ, Wilhelmy J, et al. CD8+ T cells specific for conserved coronavirus epitopes correlate with milder disease in patients with COVID-19. Science Immunology. 2021;6(61):eabg5669. doi: 10.1126/sciimmunol.abg5669.

40\. Diambra L, Alonso AM, Sookoian S, Pirola CJ. Single cell gene expression profiling of nasal ciliated cells reveals distinctive biological processes related to epigenetic mechanisms in patients with severe COVID-19. Comput Biol Med. 2022;148:105895. Epub 20220730. doi: 10.1016/j.compbiomed.2022.105895. PubMed PMID: 35926268; PubMed Central PMCID: PMCPMC9338837.

41\. Singh KK, Chaubey G, Chen JY, Suravajhala P. Decoding SARS-CoV-2 hijacking of host mitochondria in COVID-19 pathogenesis. Am J Physiol Cell Physiol. 2020;319(2):C258-c67. Epub 20200608. doi: 10.1152/ajpcell.00224.2020. PubMed PMID: 32510973; PubMed Central PMCID: PMCPMC7381712.

42\. Valdés-Aguayo JJ, Garza-Veloz I, Vargas-Rodríguez JR, Martinez-Vazquez MC, Avila-Carrasco L, Bernal-Silva S, et al. Peripheral Blood Mitochondrial DNA Levels Were Modulated by SARS-CoV-2 Infection Severity and Its Lessening Was Associated With Mortality Among Hospitalized Patients With COVID-19. Front Cell Infect Microbiol. 2021;11:754708. Epub 20211216. doi: 10.3389/fcimb.2021.754708. PubMed PMID: 34976854; PubMed Central PMCID: PMCPMC8716733.

43\. Guarnieri JW, Dybas JM, Fazelinia H, Kim MS, Frere J, Zhang Y, et al. Targeted Down Regulation Of Core Mitochondrial Genes During SARS-CoV-2 Infection. bioRxiv. 2022. Epub 20220222. doi: 10.1101/2022.02.19.481089. PubMed PMID: 35233572; PubMed Central PMCID: PMCPMC8887073.

44\. Jin X, Zhou W, Luo M, Wang P, Xu Z, Ma K, et al. Global characterization of B cell receptor repertoire in COVID-19 patients by single-cell V(D)J sequencing. Brief Bioinform. 2021;22(6). doi: 10.1093/bib/bbab192. PubMed PMID: 34015809; PubMed Central PMCID: PMCPMC8194558.

45\. Zhang Y, Yan Q, Luo K, He P, Hou R, Zhao X, et al. Analysis of B Cell Receptor Repertoires Reveals Key Signatures of the Systemic B Cell Response after SARS-CoV-2 Infection. J Virol. 2022;96(4):e0160021. Epub 20211208. doi: 10.1128/jvi.01600-21. PubMed PMID: 34878902; PubMed Central PMCID: PMCPMC8865482.

46\. Zhou JQ, Kleinstein SH. Cutting Edge: Ig H Chains Are Sufficient to Determine Most B Cell Clonal Relationships. J Immunol. 2019;203(7):1687-92. doi: 10.4049/jimmunol.1900666. PubMed PMID: 31484734.

47\. Wen W, Su W, Tang H, Le W, Zhang X, Zheng Y, et al. Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing. Cell Discov. 2020;6:31. Epub 20200504. doi: 10.1038/s41421-020-0168-9. PubMed PMID: 32377375; PubMed Central PMCID: PMCPMC7197635.

48\. Xiang H, Zhao Y, Li X, Liu P, Wang L, Wang M, et al. Landscapes and dynamic diversifications of B-cell receptor repertoires in COVID-19 patients. Human Immunology. 2022;83(2):119-29. doi: <https://doi.org/10.1016/j.humimm.2021.10.007>.

49\. Bomben R, Dal-Bo M, Benedetti D, Capello D, Forconi F, Marconi D, et al. Expression of mutated IGHV3-23 genes in chronic lymphocytic leukemia identifies a disease subset with peculiar clinical and biological features. Clin Cancer Res. 2010;16(2):620-8. Epub 20100112. doi: 10.1158/1078-0432.Ccr-09-1638. PubMed PMID: 20068100.

50\. Dal-Bo M, Del Giudice I, Bomben R, Capello D, Bertoni F, Forconi F, et al. B-cell receptor, clinical course and prognosis in chronic lymphocytic leukaemia: the growing saga of the IGHV3 subgroup gene usage. Br J Haematol. 2011;153(1):3-14. Epub 20110208. doi: 10.1111/j.1365-2141.2010.08440.x. PubMed PMID: 21303354.

51\. Pascual V, Victor K, Lelsz D, Spellerberg MB, Hamblin TJ, Thompson KM, et al. Nucleotide sequence analysis of the V regions of two IgM cold agglutinins. Evidence that the VH4-21 gene segment is responsible for the major cross-reactive idiotype. J Immunol. 1991;146(12):4385-91. PubMed PMID: 1710250.

52\. Galson JD, Schaetzle S, Bashford-Rogers RJM, Raybould MIJ, Kovaltsuk A, Kilpatrick GJ, et al. Deep Sequencing of B Cell Receptor Repertoires From COVID-19 Patients Reveals Strong Convergent Immune Signatures. Front Immunol. 2020;11:605170. Epub 20201215. doi: 10.3389/fimmu.2020.605170. PubMed PMID: 33384691; PubMed Central PMCID: PMCPMC7769841.

53\. Brouwer PJM, Caniels TG, van der Straten K, Snitselaar JL, Aldon Y, Bangaru S, et al. Potent neutralizing antibodies from COVID-19 patients define multiple targets of vulnerability. Science. 2020;369(6504):643-50. Epub 20200615. doi: 10.1126/science.abc5902. PubMed PMID: 32540902; PubMed Central PMCID: PMCPMC7299281.

54\. Tan J, Sack BK, Oyen D, Zenklusen I, Piccoli L, Barbieri S, et al. A public antibody lineage that potently inhibits malaria infection through dual binding to the circumsporozoite protein. Nat Med. 2018;24(4):401-7. Epub 20180319. doi: 10.1038/nm.4513. PubMed PMID: 29554084; PubMed Central PMCID: PMCPMC5893353.

55\. Godoy-Lozano EE, Téllez-Sosa J, Sánchez-González G, Sámano-Sánchez H, Aguilar-Salgado A, Salinas-Rodríguez A, et al. Lower IgG somatic hypermutation rates during acute dengue virus infection is compatible with a germinal center-independent B cell response. Genome Med. 2016;8(1):23. Epub 20160225. doi: 10.1186/s13073-016-0276-1. PubMed PMID: 26917418; PubMed Central PMCID: PMCPMC4766701.

56\. Ghasemi M, Kalantari S, Zubarev RA, Nafar M, Saei AA, Heidari SS, et al. Predictive Biomarker Panel in Proliferative Lupus Nephritis- Two-Dimensional Shotgun Proteomics. Iran J Kidney Dis. 2021;1(2):121-33. PubMed PMID: 33764323.

57\. Yingjuan W, Li Z, Wei C, Xiaoyuan W. Identification of prognostic genes and construction of a novel gene signature in the skin melanoma based on the tumor microenvironment. Medicine (Baltimore). 2021;100(21):e26017. doi: 10.1097/md.0000000000026017. PubMed PMID: 34032721; PubMed Central PMCID: PMCPMC8154473.

58\. Zhang JA, Zhou XY, Huang D, Luan C, Gu H, Ju M, et al. Development of an Immune-Related Gene Signature for Prognosis in Melanoma. Front Oncol. 2020;10:602555. Epub 20210121. doi: 10.3389/fonc.2020.602555. PubMed PMID: 33585219; PubMed Central PMCID: PMCPMC7874014.

59\. Prchal M, Pilz A, Simma O, Lingnau K, von Gabain A, Strobl B, et al. Type I interferons as mediators of immune adjuvants for T- and B cell-dependent acquired immunity. Vaccine. 2009;27 Suppl 6:G17-20. doi: 10.1016/j.vaccine.2009.10.016. PubMed PMID: 20006134.

60\. Cha L, Berry CM, Nolan D, Castley A, Fernandez S, French MA. Interferon-alpha, immune activation and immune dysfunction in treated HIV infection. Clin Transl Immunology. 2014;3(2):e10. Epub 20140228. doi: 10.1038/cti.2014.1. PubMed PMID: 25505958; PubMed Central PMCID: PMCPMC4232062.

61\. Boehm U, Klamp T, Groot M, Howard JC. Cellular responses to interferon-gamma. Annu Rev Immunol. 1997;15:749-95. doi: 10.1146/annurev.immunol.15.1.749. PubMed PMID: 9143706.

62\. Finkelman FD, Katona IM, Mosmann TR, Coffman RL. IFN-gamma regulates the isotypes of Ig secreted during in vivo humoral immune responses. J Immunol. 1988;140(4):1022-7. PubMed PMID: 3125247.

63\. Liu T, Zhang L, Joo D, Sun S-C. NF-κB signaling in inflammation. Signal Transduction and Targeted Therapy. 2017;2(1):17023. doi: 10.1038/sigtrans.2017.23.

64\. Zhang H, Sun SC. NF-κB in inflammation and renal diseases. Cell Biosci. 2015;5:63. Epub 20151116. doi: 10.1186/s13578-015-0056-4. PubMed PMID: 26579219; PubMed Central PMCID: PMCPMC4647710.

65\. Li F, Wang Y, Zeller KI, Potter JJ, Wonsey DR, O'Donnell KA, et al. Myc stimulates nuclearly encoded mitochondrial genes and mitochondrial biogenesis. Mol Cell Biol. 2005;25(14):6225-34. doi: 10.1128/mcb.25.14.6225-6234.2005. PubMed PMID: 15988031; PubMed Central PMCID: PMCPMC1168798.

66\. Morrish F, Hockenbery D. Myc's mastery of mitochondrial mischief. Cell Cycle. 2003;2(1):11-3. doi: 10.4161/cc.2.1.275. PubMed PMID: 12695675.

67\. Sena Laura A, Li S, Jairaman A, Prakriya M, Ezponda T, Hildeman David A, et al. Mitochondria Are Required for Antigen-Specific T Cell Activation through Reactive Oxygen Species Signaling. Immunity. 2013;38(2):225-36. doi: <https://doi.org/10.1016/j.immuni.2012.10.020>.

68\. Weinberg Samuel E, Sena Laura A, Chandel Navdeep S. Mitochondria in the Regulation of Innate and Adaptive Immunity. Immunity. 2015;42(3):406-17. doi: <https://doi.org/10.1016/j.immuni.2015.02.002>.

69\. Shenoy S. Coronavirus (Covid-19) sepsis: revisiting mitochondrial dysfunction in pathogenesis, aging, inflammation, and mortality. Inflamm Res. 2020;69(11):1077-85. Epub 20200807. doi: 10.1007/s00011-020-01389-z. PubMed PMID: 32767095; PubMed Central PMCID: PMCPMC7410962.

70\. Burtscher J, Cappellano G, Omori A, Koshiba T, Millet GP. Mitochondria: In the Cross Fire of SARS-CoV-2 and Immunity. iScience. 2020;23(10):101631. Epub 20200929. doi: 10.1016/j.isci.2020.101631. PubMed PMID: 33015593; PubMed Central PMCID: PMCPMC7524535.

71\. Aneva IY, Habtemariam S, Banach M, Sil PC, Sarkar K, Sahebnasagh A, et al. Can We Use mTOR Inhibitors for COVID-19 Therapy? Comb Chem High Throughput Screen. 2021. Epub 20211130. doi: 10.2174/1386207325666211130140923. PubMed PMID: 34847842.

72\. Husain A, Byrareddy SN. Rapamycin as a potential repurpose drug candidate for the treatment of COVID-19. Chem Biol Interact. 2020;331:109282. Epub 20201006. doi: 10.1016/j.cbi.2020.109282. PubMed PMID: 33031791; PubMed Central PMCID: PMCPMC7536130.

73\. Mullen PJ, Garcia G, Purkayastha A, Matulionis N, Schmid EW, Momcilovic M, et al. SARS-CoV-2 infection rewires host cell metabolism and is potentially susceptible to mTORC1 inhibition. Nature Communications. 2021;12(1):1876. doi: 10.1038/s41467-021-22166-4.

74\. Ramaiah MJ. mTOR inhibition and p53 activation, microRNAs: The possible therapy against pandemic COVID-19. Gene Rep. 2020;20:100765. Epub 20200630. doi: 10.1016/j.genrep.2020.100765. PubMed PMID: 32835132; PubMed Central PMCID: PMCPMC7324924.

75\. Marei HE, Althani A, Afifi N, Hasan A, Caceci T, Pozzoli G, et al. p53 signaling in cancer progression and therapy. Cancer Cell International. 2021;21(1):703. doi: 10.1186/s12935-021-02396-8.

76\. Cardozo CM, Hainaut P. Viral strategies for circumventing p53: the case of severe acute respiratory syndrome coronavirus. Curr Opin Oncol. 2021;33(2):149-58. doi: 10.1097/cco.0000000000000713. PubMed PMID: 33405482; PubMed Central PMCID: PMCPMC7924916.

77\. Ma-Lauer Y, Carbajo-Lozoya J, Hein MY, Müller MA, Deng W, Lei J, et al. p53 down-regulates SARS coronavirus replication and is targeted by the SARS-unique domain and PL\<sup\>pro\</sup\> via E3 ubiquitin ligase RCHY1. Proceedings of the National Academy of Sciences. 2016;113(35):E5192-E201. doi: doi:10.1073/pnas.1603435113.

78\. Gralinski LE, Sheahan TP, Morrison TE, Menachery VD, Jensen K, Leist SR, et al. Complement Activation Contributes to Severe Acute Respiratory Syndrome Coronavirus Pathogenesis. mBio. 2018;9(5). Epub 20181009. doi: 10.1128/mBio.01753-18. PubMed PMID: 30301856; PubMed Central PMCID: PMCPMC6178621.

## 
