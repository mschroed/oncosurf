library(Biobase)
library(affy)                    

#load("data/minn2007.RData")
#load("data/nki2002.RData")
#load("data/sotiriou2006.RData")
#load("data/miller2005.RData")
#load("data/transbig2006affy.RData")
#load("data/schmidt2008.RData")   
##> intersect(colnames(pData(vdx)),intersect(colnames(pData(upp)),intersect(colnames(pData(unt)),intersect(colnames(pData(transbig)),intersect(colnames(pData(mainz)),colnames(pData(nki)))))))
## [1] "samplename"    "series"        "dataset"       "filename"     
## [5] "id"            "e.dmfs"        "t.dmfs"        "node"         
## [9] "t.rfs"         "e.rfs"         "er"            "size"         
##[13] "age"           "grade"         "treatment"     "tissue"       
##[17] "t.os"          "e.os"          "pgr"           "her2"         
##[21] "brca.mutation"

## set the rank of columnnames with the names from the intersection of all 6 eSets
columnRank <- c("samplename","dataset","series","id","filename","size","age","er","grade","pgr","her2","brca.mutation","e.dmfs","t.dmfs","node","t.rfs","e.rfs","treatment","tissue","t.os","e.os")

## creation of the first eSet
load("data/minn2007.RData")
demo$T <- demo$size
demo <- demo[,!is.element(colnames(demo), "size.cat")]
demo <- demo[,!is.element(colnames(demo), "primary.bc.unt")]
demo$size <- c(rep(NA,length(demo$size)))
demo <- demo[,columnRank]
#annt=phenotyp daten
metadata<-data.frame(labelDescription=colnames(demo), row.names=colnames(demo))
phenoD<-new("AnnotatedDataFrame", data=demo, varMetadata=metadata)
## probe annotations
metadata<-data.frame(labelDescription=colnames(annot), row.names=colnames(annot))
featureD <- new("AnnotatedDataFrame", data=annot, varMetadata=metadata)
experimentD <- new("MIAME",
  name = "VDX",
  lab = "Veridex LLC, a Johnson & Johnson Company, San Diego, CA, USA. Department of Radiaton and Cellular Oncology, Center for Molecular Oncology, and Ludwig Center for Metastasis Research, University of Chicago, Chicago, IL, USA",
  contact = "Dr John Foekens, Erasmus MC Josephine Nefkens Institute, Netherlands: <j.foekens@erasmusmc.nl>, Joan Massague,  <j-massagie@ski.mskcc.org>",
  title = "Gene-expression profiles to predict distant metastasis of kymph-node-negative primary breast cancer. Lung metastasis genes coule breast tumor size and metastatic spread.",
  abstract = "Wang et al. 2005 and Minn et al. 2007. The association between large tumor size and metastatic risk in a majority of clinical cancers has led to questions as to whether these observations are causally related or whether one is simply a marker for the other. This is partly due to an uncertainty about how metastasis-promoting gene expression changes can arise in primary tumors. We investigated this question through the analysis of a previously defined 'lung metastasis gene-expression signa- ture'(LMS) that mediates experimental breast cancer metastasis selectively to the lung and is expressed by primary human breast cancer with a high risk for developing lung metastasis. Experimentally, we demonstrate that the LMS promotes primary tumor growth that enriches for LMS cells, and it allows for intravasation after reaching a critical tumor size. Clinically, this corresponds to LMS tumors being larger at diagnosis compared with LMS tumors and to a marked rise in the incidence of metastasis after LMS tumors reach 2 cm. Patients with LMS-expressing primary tumors selectively fail in the lung compared with the bone or other visceral sites and have a worse overall survival. The mechanistic linkage between metastasis gene expression, accelerated tumor growth, and likelihood of metastatic recurrence provided by the LMS may help to explain observations of prognostic gene signatures in primary cancer and how tumor growth can both lead to metastasis and be a marker for cells destined to metastasize.",
  url = "GEO accession number: GSE2034 & GSE5327 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034>, <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5327>",
  pubMedIds = "17420468")
  
arrayT <- "hgu133a"
#arraytype is an string
vdx <- new("ExpressionSet",
          exprs=t(data),
          phenoData=phenoD,
          featureData=featureD,          
          annotation=arrayT,
          experimentData=experimentD)



## creation of the second eSet
load("data/nki2002.RData")
demo <- demo[,!is.element(colnames(demo), "gg")]
demo <- demo[,columnRank]
#annt=phenotyp daten
metadata<-data.frame(labelDescription=colnames(demo), row.names=colnames(demo))
phenoD<-new("AnnotatedDataFrame", data=demo, varMetadata=metadata)
## probe annotations
metadata<-data.frame(labelDescription=colnames(annot), row.names=colnames(annot))
featureD <- new("AnnotatedDataFrame", data=annot, varMetadata=metadata)
experimentD <- new("MIAME",
  name = "NKI",
  lab = "Divisions of Diagnostic Oncology, Radiotherapy and Molecular Carvinogenesis and Center for Biomedical Genetics, The Netherland Cancer Institute, Amsterdam, The Netherlands.",
  contact = "Stephen H. Friend <stephen_friend@merck.com>",
  title = "Gene expression profiling predicts clinical outcome of breast cancer. A gene-expresion signature as a predictor  of survival in breast cancer.",
  abstract = "Van de Vijver et al. 2002 and Laura J. van't Veer 2002. Series of 295 concescutive patients with primary breast carcinomas as having a gene-expression signature associated with either a poor prognosis or a good prognosis. All patients had stage I or II breast cancer and were younger than 53 years old; 151 had lymph-node-negaitve disease, and 144 had lymph-node-positive disease.",
  url = "http://www.rii.com/publications/2002/vantveer.html")
arrayT <- "rosetta"
#arraytype is an string
nki <- new("ExpressionSet",
          exprs=t(data),
          phenoData=phenoD,
          featureData=featureD,
          annotation=arrayT,
          experimentData=experimentD)


## creation of the third eSet
load("data/sotiriou2006.RData")
demo <- demo[,!is.element(colnames(demo), "gg")]
demo <- demo[,columnRank]
#annt=phenotyp daten
metadata<-data.frame(labelDescription=colnames(demo), row.names=colnames(demo))
phenoD<-new("AnnotatedDataFrame", data=demo, varMetadata=metadata)
## probe annotations
metadata<-data.frame(labelDescription=colnames(annot), row.names=colnames(annot))
featureD <- new("AnnotatedDataFrame", data=annot, varMetadata=metadata)
experimentD <- new("MIAME",
  name = "UNT",
  lab = "Functional Genomics and Translational Research Unit, Jules Bordet Institute, Universite Libre de Bruxelles, Brussels, Belgium.",
  contact = "Christos Sotiriou <christos.sotiriou@bordet.be>",
  title = "Gene expression profiling in breast cancer: understanding the molecular basis of histologic grade to improve prognosis.",
  abstract = "Sotiriou et al. 2006. Background: Histologic grade in breast cancer provides clinically important prognostic information. However, 30%-60% of tumors are classified as histologic grade 2. This grade is associated with an intermediate risk of recurrence and is thus not informative for clinical decision making. We examined whether histologic grade was associated with gene expression profiles of breast cancers and whether such profiles could be used to improve histologic grading. Methods: We analyzed microarray data from 189 invasive breast carcinomas and from three published gene expression datasets from breast carcinomas. We identified differentially expressed genes in a training set of 64 estrogen receptor (ER)-positive tumor samples by comparing expression profiles between histologic grade 3 tumors and histologic grade 1 tumors and used the expression of these genes to define the gene expression grade index. Data from 597 independent tumors were used to evaluate the association between relapse-free survival and the gene expression grade index in a Kaplan-Meier analysis. All statistical tests were two-sided. Results: We identified 97 genes in our training set that were associated with histologic grade; most of these genes were involved in cell cycle regulation and proliferation. In validation datasets, the gene expression grade index was strongly associated with histologic grade 1 and 3 status; however, among histologic grade 2 tumors, the index spanned the values for histologic grade 1-3 tumors. Among patients with histologic grade 2 tumors, a high gene expression grade index was associated with a higher risk of recurrence than a low gene expression grade index (hazard ratio = 3.61, 95% confidence interval = 2.25 to 5.78; P < .001, log-rank test). Conclusions: Gene expression grade index appeared to reclassify patients with histologic grade 2 tumors into two groups with high versus low risks of recurrence. This approach may improve the accuracy of tumor grading and thus its prognostic value.",
  pubMedIds = "16478745",
  url = "GEO accession number: GEO GSE2990 & GSE6532 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2990>, <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532>")
arrayT <- "hgu133ab"
#arraytype is an string
unt <- new("ExpressionSet",
          exprs=t(data),
          phenoData=phenoD,
          featureData=featureD,
          annotation=arrayT,
          experimentData=experimentD)



## creation of the forth eSet
load("data/miller2005.RData")
demo <- demo[,!is.element(colnames(demo), "gg")]
demo <- demo[,!is.element(colnames(demo), "treatment2")]
demo <- demo[,columnRank]
#annt=phenotyp daten
metadata<-data.frame(labelDescription=colnames(demo), row.names=colnames(demo))
phenoD<-new("AnnotatedDataFrame", data=demo, varMetadata=metadata)
## probe annotations
metadata<-data.frame(labelDescription=colnames(annot), row.names=colnames(annot))
featureD <- new("AnnotatedDataFrame", data=annot, varMetadata=metadata)
experimentD <- new("MIAME",
  name = "UPP",
  lab = "Genome Institute of Singapore, Singapore and  Department of Oncology and Pathology, Radiumhemmet, Karolinksa Institue and Hospital, Stockholm, Sweden",
  contact = "Lance D. Miller <millerl@gis.a-star.edu.sg> and Edison T. Liu <luie@gis.a-star.edu.sg>",
  title = "An expression signature for p53 status in human breast caner predicts mutation status, transcriptional effects, and patient survival.",
  abstract = "Miller et al. 2005. Perturbations of the p53 pathway are associated with more aggressive and therapeutically refractory tumors. However, molecular assessment of p53 status, by using sequence analysis and immunohistochemistry, are incomplete assessors of p53 functional effects. We posited that the transcriptional fingerprint is a more definitive downstream indicator of p53 function. Herein, we analyzed transcript profiles of 251 p53-sequenced primary breast tumors and identified a clinically embedded 32-gene expression signature that distinguishes p53-mutant and wild-type tumors of different histologies and outperforms sequence-based assessments of p53 in predicting prognosis and therapeutic response. Moreover, the p53 signature identified a subset of aggressive tumors absent of sequence mutations in p53 yet exhibiting expression characteristics consistent with p53 deficiency because of attenuated p53 transcript levels. Our results show the primary importance of p53 functional status in predicting clinical breast cancer behavior.",
  url = "GEO accession number: GSE3494 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3494>",
  pubMedIds = "16141321") 
arrayT <- "hgu133ab"
## create eSet
upp <- new("ExpressionSet",
          exprs=t(data),
          phenoData=phenoD,
          featureData=featureD,
          annotation=arrayT,
          experimentData=experimentD)
          

## creation of the fifth eSet
load("data/transbig2006affy.RData")
demo <- demo[,!is.element(colnames(demo), "gg")]
demo <- demo[,!is.element(colnames(demo), "veridex_score")]
demo <- demo[,!is.element(colnames(demo), "agendia_score")]
demo <- demo[,!is.element(colnames(demo), "Date_of_Birth")]
demo <- demo[,!is.element(colnames(demo), "Date_of_Diagnosis")]
demo <- demo[,columnRank]
#annt=phenotyp daten
metadata<-data.frame(labelDescription=colnames(demo), row.names=colnames(demo))
phenoD<-new("AnnotatedDataFrame", data=demo, varMetadata=metadata)
## probe annotations
metadata<-data.frame(labelDescription=colnames(annot), row.names=colnames(annot))
featureD <- new("AnnotatedDataFrame", data=annot, varMetadata=metadata)
experimentD <- new("MIAME",
  name = "TRANSBIG",
  lab = "Institute Jules Bordet, Universite Libre de Bruxelles, Brussels, Belgium",
  contact = "TRANSBIG Consortium <transbig@bordet.be>",
  title = "Strong Time dependence of the 76-gene prognostic signature for node-negative breast cancer patients int the TRANSBIG multicenter independent validation series",
  abstract = "Desmedt et al. 2007. Recently a 76-gene prognostic signature able to predict distant metastases in lymph node-negative (N-) breast cancer patients was reported. The aims of this study conducted by TRANSBIG were to independently validate these results and to compare the outcome with clinical risk assessment. Materials and Methods: Gene expression profiling of frozen samples from 198 N- systemically untreated patients was performed at the Bordet Institute, blinded to clinical data and independent of Veridex. Genomic risk was defined by Veridex, blinded to clinical data. Survival analyses, done by an independent statistician, were performed with the genomic risk and adjusted for the clinical risk, defined by Adjuvant!Online. Results: The actual 5- and 10-year time to distant metastasis (TDM) were 98% (88%-100%) and 94% (83%-98%) respectively for the good profile group and 76% (68%- 82%) and 73% (65%-79%) for the poor profile group. The actual 5- and 10-year overall survival (OS) were 98% (88%-100%) and 87% (73%-94%) respectively for the good profile group and 84% (77%-89%) and 72% (63%-78%) for the poor profile group. We observed a strong time-dependency of this signature, leading to an adjusted HR of 13.58 (1.85-99.63) and 8.20 (1.10-60.90) at 5 years, and 5.11 (1.57-16.67) and 2.55 (1.07-6.10) at 10 years for TDM and OS respectively. Conclusion: This independent validation confirmed the performance of the 76-gene signature and adds to the growing evidence that gene expression signatures are of clinical relevance, especially for identifying patients at high risk of early distant metastases.",
  url = "GEO accession number: GSE7390 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7390>",
  pubMedIds = "17545524")
arrayT <- "hgu133a"  
## create eSet
transbig <- new("ExpressionSet",
          exprs=t(data),
          phenoData=phenoD,
          featureData=featureD,          
          annotation=arrayT,
          experimentData=experimentD)
          
          
## creation of the sixth eSet
load("data/schmidt2008.RData")
demo <- demo[,columnRank]
#annt=phenotyp daten
metadata<-data.frame(labelDescription=colnames(demo), row.names=colnames(demo))
phenoD<-new("AnnotatedDataFrame", data=demo, varMetadata=metadata)
## probe annotations
metadata<-data.frame(labelDescription=colnames(annot), row.names=colnames(annot))
featureD <- new("AnnotatedDataFrame", data=annot, varMetadata=metadata)
experimentD <- new("MIAME",
  name = "MAINZ",
  lab = "Department of Obstetrics and Gynecology, Medical School, Johannes Gutenberg University, Mainz, Germany",
  contact = "Mathias Gehrmann <mathias.gehrmann@siemens.com>",
  title = "The humoral immune system has a key prognostic impact in node-negative breast cancer.",
  abstract = "Schmidt et al. 2008. Background: Estrogen receptor (ER) expression and proliferative activity are established prognostic factors in breast cancer. In a search for additional prognostic motifs, we analyzed the gene expression patterns of 200 tumors of patients who were not treated by systemic therapy after surgery using a discovery approach. After performing hierarchical cluster analysis, we identified coregulated genes related to the biological process of proliferation, steroid hormone receptor expression, as well as B-cell and T-cell infiltration. We calculated metagenes as a surrogate for all genes contained within a particular cluster and visualized the relative expression in relation to time to metastasis with principal component analysis. Distinct pat- terns led to the hypothesis of a prognostic role of the immune system in tumors with high expression of proliferation- associated genes. In multivariate Cox regression analysis, the proliferation metagene showed a significant association with metastasis-free survival of the whole discovery cohort [hazard ratio (HR), 2.20; 95% confidence interval (95% CI), 1.40-3.46]. The B-cell metagene showed additional indepen- dent prognostic information in carcinomas with high prolif- erative activity (HR, 0.66; 95% CI, 0.46-0.97). A prognostic influence of the B-cell metagene was independently confirmed by multivariate analysis in a first validation cohort enriched for high-grade tumors (n = 286; HR, 0.78; 95% CI, 0.62-0.98) and a second validation cohort enriched for younger patients (n = 302; HR, 0.83; 95% CI, 0.7-0.97). Thus, we could show in three cohorts of untreated, node-negative breast cancer patients that the humoral immune system plays a pivotal role in metastasis-free survival of carcinomas of the breast.",
  url = "GEO accession number: GSE11121 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121>",
  pubMedIds = "18593943")
arrayT <- "hgu133a"    
## create eSet  
mainz <- new("ExpressionSet",
          exprs=t(data),
          phenoData=phenoD,
          featureData=featureD,          
          annotation=arrayT,
          experimentData=experimentD)
          
##package.skeleton("datasets")
