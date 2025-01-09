#' @rdname Datasets
#' @name Datasets
#' @title Datasets
#' @description The datasets are based on the official .zip file. A table for the dataset names and file names sorted by file name is here:
#' \tabular{ll}{
#' \bold{Dataset name} \tab \bold{File name} \cr
#' -------------------- \tab ----------------- \cr
#' \link{illustration} \tab "A numerical illustration.dat" \cr
#' \link{leukaemia} \tab "Bone marrow transplantation in the treatment of leukaemia.dat" \cr
#' \link{bone_marrow} \tab "Bone marrow transplantation.dat" \cr
#' \link{ovarian} \tab "Chemotherapy in ovarian cancer patients.dat" \cr
#' \link{active_hepatitis} \tab "Chronic active hepatitis.dat" \cr
#' \link{granulomatous} \tab "Chronic granulomatous disease.dat" \cr
#' \link{tamoxifen} \tab "Clinical trial of tamoxifen in breast cancer patients.dat" \cr
#' \link{prostatic} \tab "Comparison of two treatments for prostatic cancer.dat" \cr
#' \link{kidneytx} \tab "Comparisons between kidney transplant centres.dat" \cr
#' \link{liverbase} \tab "Data from a cirrhosis study (baseline).dat" \cr
#' \link{liver_counting} \tab "Data from a cirrhosis study (in counting process format).dat" \cr
#' \link{lbrdata0} \tab "Data from a cirrhosis study (lbr data).dat" \cr
#' \link{HELP} \tab "Health evaluation and linkage to primary care.dat" \cr
#' \link{dialysis} \tab "Infection in patients on dialysis.dat" \cr
#' \link{bone_marrow_tx} \tab "Patient outcome following bone marrow transplantation.dat" \cr
#' \link{bcancer} \tab "Prognosis for women with breast cancer.dat" \cr
#' \link{pulmonary} \tab "Pulmonary metastasis.dat" \cr
#' \link{breast_rfs} \tab "Recurrence free survival in breast cancer patients.dat" \cr
#' \link{ulcer} \tab "Recurrence of an ulcer.dat" \cr
#' \link{bladder} \tab "Recurrence of bladder cancer.dat" \cr
#' \link{mammary} \tab "Recurrence of mammary tumours in female rats.dat" \cr
#' \link{valve} \tab "Survival following aortic valve replacement.dat" \cr
#' \link{tplant} \tab "Survival following kidney transplantation.dat" \cr
#' \link{ducks} \tab "Survival of black ducks.dat" \cr
#' \link{mice} \tab "Survival of laboratory mice.dat" \cr
#' \link{liver} \tab "Survival of liver transplant recipients.dat" \cr
#' \link{myeloma} \tab "Survival of multiple myeloma patients.dat" \cr
#' \link{lung} \tab "Survival of patients registered for a lung transplant.dat" \cr
#' \link{gcancer} \tab "Survival of patients with gastric cancer.dat" \cr
#' \link{melanoma} \tab "Survival times of patients with melanoma .dat" \cr
#' \link{livertx} \tab "Time to death while waiting for a liver transplant.dat" \cr
#' \link{IUD} \tab "Time to discontinuation of the use of an IUD.dat" \cr
#' \link{kidney} \tab "Treatment of hypernephroma.dat" \cr
#' }
#' And now sorted by the dataset names:
#' \tabular{ll}{
#' \bold{Dataset name} \tab \bold{File name} \cr
#' -------------------- \tab ----------------- \cr
#' \link{active_hepatitis} \tab "Chronic active hepatitis.dat" \cr
#' \link{bcancer} \tab "Prognosis for women with breast cancer.dat" \cr
#' \link{bladder} \tab "Recurrence of bladder cancer.dat" \cr
#' \link{bone_marrow} \tab "Bone marrow transplantation.dat" \cr
#' \link{bone_marrow_tx} \tab "Patient outcome following bone marrow transplantation.dat" \cr
#' \link{breast_rfs} \tab "Recurrence free survival in breast cancer patients.dat" \cr
#' \link{dialysis} \tab "Infection in patients on dialysis.dat" \cr
#' \link{ducks} \tab "Survival of black ducks.dat" \cr
#' \link{gcancer} \tab "Survival of patients with gastric cancer.dat" \cr
#' \link{granulomatous} \tab "Chronic granulomatous disease.dat" \cr
#' \link{HELP} \tab "Health evaluation and linkage to primary care.dat" \cr
#' \link{illustration} \tab "A numerical illustration.dat" \cr
#' \link{IUD} \tab "Time to discontinuation of the use of an IUD.dat" \cr
#' \link{kidney} \tab "Treatment of hypernephroma.dat" \cr
#' \link{kidneytx} \tab "Comparisons between kidney transplant centres.dat" \cr
#' \link{lbrdata0} \tab "Data from a cirrhosis study (lbr data).dat" \cr
#' \link{leukaemia} \tab "Bone marrow transplantation in the treatment of leukaemia.dat" \cr
#' \link{liver} \tab "Survival of liver transplant recipients.dat" \cr
#' \link{liver_counting} \tab "Data from a cirrhosis study (in counting process format).dat" \cr
#' \link{liverbase} \tab "Data from a cirrhosis study (baseline).dat" \cr
#' \link{livertx} \tab "Time to death while waiting for a liver transplant.dat" \cr
#' \link{lung} \tab "Survival of patients registered for a lung transplant.dat" \cr
#' \link{mammary} \tab "Recurrence of mammary tumours in female rats.dat" \cr
#' \link{melanoma} \tab "Survival times of patients with melanoma .dat" \cr
#' \link{mice} \tab "Survival of laboratory mice.dat" \cr
#' \link{myeloma} \tab "Survival of multiple myeloma patients.dat" \cr
#' \link{ovarian} \tab "Chemotherapy in ovarian cancer patients.dat" \cr
#' \link{prostatic} \tab "Comparison of two treatments for prostatic cancer.dat" \cr
#' \link{pulmonary} \tab "Pulmonary metastasis.dat" \cr
#' \link{tamoxifen} \tab "Clinical trial of tamoxifen in breast cancer patients.dat" \cr
#' \link{tplant} \tab "Survival following kidney transplantation.dat" \cr
#' \link{ulcer} \tab "Recurrence of an ulcer.dat" \cr
#' \link{valve} \tab "Survival following aortic valve replacement.dat" \cr
#' }
#' @source \url{https://s3-eu-west-1.amazonaws.com/s3-euw1-ap-pe-ws4-cws-documents.ri-prod/9781032252858/Data%20sets%20from%20Modelling%20Survival%20Data%20in%20Medical%20Research%2C%204th%20edition.zip}
"_PACKAGE"

#' @title A numerical illustration
#' @description Artificial data on patient survival classified according to factors a and b
#' @format A data frame with 37 rows and 4 variables:
#' \describe{
#'   \item{\code{a}}{integer factor a}
#'   \item{\code{b}}{integer factor b}
#'   \item{\code{time}}{integer event time}
#'   \item{\code{status}}{integer event status (1=event, 0=right censored)} 
#'}
#' @details See Collett (2023). 
"illustration"

#' @title Bone marrow transplantation in the treatment of leukaemia
#' @description Bone marrow transplantation in the treatment of leukaemia
#' @format A data frame with 23 rows and 8 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{time}}{integer survival time in days}
#'   \item{\code{status}}{integer event indicator (0=alive, 1=dead)}
#'   \item{\code{group}}{integer disease group (1=ALL, 2=low-risk AML, 3=high-risk AML)}
#'   \item{\code{page}}{integer age of patient in years}
#'   \item{\code{dage}}{integer age of donor in years}
#'   \item{\code{precovery}}{integer platelet recovery indicator (0=no, 1=yes)}
#'   \item{\code{ptime}}{character time in days to return of platelets to normal level (if precovery=1)} 
#'}
#' @details See Collett (2023). Note that ptime will need conversion:).
"leukaemia"

#' @title Bone marrow transplantation
#' @description A study of 37 patients with leukaemia in complete remission who received a non-depleted allogenic bone marrow transplant.
#' @format A data frame with 37 rows and 9 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient number (1-37)}
#'   \item{\code{time}}{integer survival time in days}
#'   \item{\code{status}}{integer status of patient (0=alive, 1=dead)}
#'   \item{\code{rage}}{integer age of patient in years}
#'   \item{\code{dage}}{integer age of donor in years}
#'   \item{\code{type}}{integer type of leukaemia (1=AML, 2=ALL, 3=CML)}
#'   \item{\code{preg}}{integer Donor pregnancy (0=no, 1=yes)}
#'   \item{\code{index}}{double index of cell-lymphocyte reactions}
#'   \item{\code{gvhd}}{integer graft-versus-host disease (0=no, 1=yes)} 
#'}
#' @details See Collett (2023)
"bone_marrow"

#' @title Chemotherapy in ovarian cancer patients
#' @description Trial for treatment of ovarian cancer patients comparing cyclophosphamide alone with cyclophosphamide combined with adriamycin.
#' @format A data frame with 26 rows and 7 variables:
#' \describe{
#'   \item{\code{patient}}{integer identifer}
#'   \item{\code{time}}{integer survival time from randomisation in days}
#'   \item{\code{status}}{integer event indicator (0=right censored, 1=event)}
#'   \item{\code{treat}}{integer treatment (1=single, 2=combined)}
#'   \item{\code{age}}{integer age of patients in years}
#'   \item{\code{rdisease}}{integer extent of residual disease (1=incomplete, 2=complete)}
#'   \item{\code{perf}}{integer performance status (1=good, 2=poor)} 
#'}
#' @details See Collett (2023)
"ovarian"

#' @title Chronic active hepatitis
#' @description Clinical trial of 44 patients with chronic active hepatitis randomised to either the drug prednisolone or an untreated control group.
#' @format A data frame with 44 rows and 3 variables:
#' \describe{
#'   \item{\code{treatment}}{integer treatment (1=prednisolone, 2=control)}
#'   \item{\code{time}}{integer survival time from admission to study (months)}
#'   \item{\code{status}}{integer event indicator (1=event, 0=right censored)} 
#'}
#' @details See Collett (2023)
"active_hepatitis"

#' @title Chronic granulomatous disease
#' @description Trial comparing interferon with a placebo.
#' @format A data frame with 128 rows and 12 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient number (1-128)}
#'   \item{\code{time}}{integer time to first infection (days)}
#'   \item{\code{status}}{integer status of patient (0=censored, 1=infection)}
#'   \item{\code{centre}}{integer treatment centre; see Collett (2023, page 504)}
#'   \item{\code{treat}}{integer treatment group (0=placebo, 1=interferon)}
#'   \item{\code{age}}{integer age in years}
#'   \item{\code{sex}}{integer sex (1=male, 2=female)}
#'   \item{\code{height}}{double height in cm}
#'   \item{\code{weight}}{double weight in kg}
#'   \item{\code{pattern}}{integer pattern of inheritance (1=X-linked, 2=autosomal recessive)}
#'   \item{\code{cort}}{integer use of corticosteroids at trial entry (1=used, 2=not used)}
#'   \item{\code{anti}}{integer Use of antibiotics at trial entry (1=used, 2=not used)} 
#'}
#' @details See Collett (2023)
"granulomatous"

#' @title Clinical trial of tamoxifen in breast cancer patients
#' @description Clinical trial for breast cancer patients comparing combined tamoxifen and radiotherapy with tamoxifen alone.
#' @format A data frame with 641 rows and 18 variables:
#' \describe{
#'   \item{\code{id}}{integer patient identifier}
#'   \item{\code{treat}}{integer treatment group (0=tamoxifen+radiotherapy, 1=tamoxifen)}
#'   \item{\code{age}}{integer patient age at study entry (years)}
#'   \item{\code{size}}{double tumour size (cm)}
#'   \item{\code{hist}}{integer tumour histology (1=ductal, 2=lobular, 3=medullary, 4=mixed, 5=other)}
#'   \item{\code{hr}}{integer hormone receptor level (0=negative, 1=positive)}
#'   \item{\code{hb}}{integer Haemoglobin level (g/l)}
#'   \item{\code{andis}}{integer axillary relapse (0=no, 1=yes)}
#'   \item{\code{lsurv}}{integer time to local relapse or last follow-up (days)}
#'   \item{\code{ls}}{integer local relapse (0=no, 1=yes))}
#'   \item{\code{asurv}}{integer time to axillary relapse or last follow-up (days)}
#'   \item{\code{as}}{integer axillary relapse (0=no, 1=yes)}
#'   \item{\code{dsurv}}{integer Time to distant relapse or last follow-up (days)}
#'   \item{\code{ds}}{integer distant relapse (0=no, 1=yes)}
#'   \item{\code{msurv}}{integer time to second malignancy or last follow-up (days)}
#'   \item{\code{ms}}{integer second malignancy (0=no, 1=yes)}
#'   \item{\code{tsurv}}{integer time from randomisation to death or last follow-up (days)}
#'   \item{\code{ts}}{integer status at last follw-up (0=alive, 1=dead)} 
#'}
#' @details See Collett (2023)
"tamoxifen"

#' @title Comparison of two treatments for prostatic cancer
#' @description Randomised controlled trial from the Veteran's Administration Cooperative Urological Research Group. Includes patients who had stage III cancers and were randomised to placebo or daily oral treatment with 1.0 mg of diethylstilbesterol (DES).
#' @format A data frame with 38 rows and 8 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient identifier}
#'   \item{\code{treatment}}{integer treatment indicator (1=placebo; 2=daily treatment with 1.0 mg of diethylstilbesterol (DES))}
#'   \item{\code{time}}{integ
#' er survival time from trial entry to end of follow-up in months}
#'   \item{\code{status}}{integer for follow-up status (0=alive or died from other causes, 1=died from prostate cancer}
#'   \item{\code{age}}{integer age at trial entry in years}
#'   \item{\code{shb}}{double serum hemoglobin at trial entry in g/dL}
#'   \item{\code{size}}{integer size of the primary tumour in cm^3}
#'   \item{\code{index}}{integer Gleason index based on histopathology} 
#'}
#' @details TBC.
#' @references
#' Andrews DF, Herzberg AM. Data: a collection of problems from many fields for the student and research worker. Springer Series in Statistics; Springer New York, NY; 1985. \doi{10.1007/978-1-4612-5098-2}
"prostatic"

#' @title Comparisons between kidney transplant centres
#' @description Transplant survival rates by recipients of organs from deceased donors. No event was defined as being alive  with a functioning graft at the last known follow-up.
#' @format A data frame with 1439 rows and 9 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{centre}}{integer transplant centre (1-8)}
#'   \item{\code{tsurv}}{integer transplant survival time (days)}
#'   \item{\code{tcens}}{integer event indicator (0=censored, 1=transplant failure)}
#'   \item{\code{dage}}{integer donor age (years)}
#'   \item{\code{dtype}}{integer donor type (0=deceased following brain death, 1=circulatory death)}
#'   \item{\code{rage}}{integer recipient age (years)}
#'   \item{\code{diab}}{integer diabetic status (0=absent, 1=present)}
#'   \item{\code{cit}}{double cold ischaemic time (hours)} 
#'}
#' @details See Collett (2023). Thirty-five patients had tsurv==0 (that is, the transplanted kidney did not function).
"kidneytx"

#' @title Data from a cirrhosis study (baseline)
#' @description Articial data
#' @format A data frame with 12 rows and 6 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{time}}{integer survival time in days}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=uncensored)}
#'   \item{\code{age}}{integer age of the patient (years)}
#'   \item{\code{treat}}{integer treatment group (0=placebo, 1=Liverol)}
#'   \item{\code{lbr}}{double logarithm of bilirubin level} 
#'}
#' @details See Collett (2023)
"liverbase"

#' @title Data from a cirrhosis study (in counting process format)
#' @description Artificial data
#' @format A data frame with 54 rows and 7 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{start}}{integer start time (days)}
#'   \item{\code{stop}}{integer stop time (days)}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=uncensored)}
#'   \item{\code{treat}}{integer treatment group (0=placebo, 1=Liverol)}
#'   \item{\code{age}}{integer age of the patient at start of study (years)}
#'   \item{\code{lbrt}}{double logarithm of bilirubin level} 
#'}
#' @details See Collett (2023). Note that the variable for log of bilirubin differs to that for "liverbase".
"liver_counting"

#' @title Data from a cirrhosis study (lbr data)
#' @description DATASET_DESCRIPTION
#' @format A data frame with 42 rows and 3 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{time}}{integer date of measurement (days)}
#'   \item{\code{lbr}}{double log bilirubin level} 
#'}
#' @details See Collett (2023)
"lbrdata0"

#' @title Health evaluation and linkage to primary care
#' @description A clinical trial for patients in a residential detoxification programme. Patients were randomised to either get a referral to a HELP clinic or not.
#' @format A data frame with 447 rows and 7 variables:
#' \describe{
#'     \item{\code{subject}}{integer subject id}
#'     \item{\code{days}}{integer time to linkage to primary care in days}
#'     \item{\code{status}}{integer event indicator (0=no linkage, 1=linkage)}
#'     \item{\code{age}}{integer age of patient in years}
#'     \item{\code{gender}}{integer gender of the patient (0=female, 1=male)}
#'     \item{\code{housing}}{integer Homelessness status (0=homeless, 1=housed)}
#'     \item{\code{linkage}}{integer assistance to linking to healthcare (0=no, 1=yes)}
#' }
#' @details Collett (2023) defines this dataset as "help", however
#'     that leads to issues with using R's help system. We have
#'     changed the dataset name to "HELP". Moreover, the book uses the variables "Time" an d"Help", whereas the dataset includes variables "days" and "linkage", respectively.
"HELP"

#' @title Infection in patients on dialysis
#' @description Time from dialysis to infection for patients with diseases of the kidney.
#' @format A data frame with 13 rows and 5 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{time}}{integer time to infection (days)}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=infection)}
#'   \item{\code{age}}{integer age in years}
#'   \item{\code{sex}}{integer sex of the patient (1=male, 2=female)} 
#'}
#' @details See Collett (2023)
"dialysis"

#' @title Patient outcome following bone marrow transplantation
#' @description Patient outcome following bone marrow transplantation
#' @format A data frame with 2204 rows and 9 variables:
#' \describe{
#'   \item{\code{id}}{integer patient id}
#'   \item{\code{leukaemia}}{character type of leukaemia (CML,ALL,AML)}
#'   \item{\code{age}}{character age group of patient in years (<=20, 21-40, >40))}
#'   \item{\code{match}}{integer indicator for whether there was a donor gender match (0=no, 1=yes)}
#'   \item{\code{tcell}}{integer indicator for whether there was T-cell depletion (1=yes, n=no)}
#'   \item{\code{ptime}}{integer time to platelet recovery (days)}
#'   \item{\code{pcens}}{integer event indicator for platelet recovery (1=event, 0=censored)}
#'   \item{\code{rdtime}}{integer time to relapse of death (days)}
#'   \item{\code{rdcens}}{integer event indicator for relapse or death (1=event, 0=censored)} 
#'}
#' @details See Collett (2023)
"bone_marrow_tx"

#' @title Prognosis for women with breast cancer
#' @description For female breast cancer patients from Middlesex
#'     Hospital. The dataset includes the result of staining using Helix pomatia agglutinin (HPA).
#' @format A data frame with 45 rows and 3 variables:
#' \describe{
#'     \item{\code{stain}}{integer for negative staining (=1) or
#'     positive staining (=2)}
#'     \item{\code{time}}{integer time in
#'     months for survival}
#'     \item{\code{status}}{integer for status at
#'     end of follow-up (0=censored, 1=death)}
#' }
#' @details
#' For details about the study design, see Leathem and Brooks (1987).
#'
#' The dataset is described in Example 1.2 and Table 1.2 (Collett, 2023, pages 6-7).
#' @examples
#' library(survival)
#' plot(survfit(Surv(time,status)~stain, data=bcancer), col=1:2, xlab="Survival time (months)",
#' ylab="Survival")
#' legend("topright", legend=c("Negative staining","Positive staining"), col=1:2, lty=1,
#' bty="n")
#' @references
#' Leathem AJ, Brooks S. Predictive value of lectin binding on breast-cancer recurrence and survival. The Lancet. 1987 May 9;329(8541):1054-6. \doi{10.1016/S0140-6736(87)90482-X}
"bcancer"

#' @title Pulmonary metastasis
#' @description A very simple dataset with no censoring
#' @format A data frame with 11 rows and 1 variables:
#' \describe{
#'   \item{\code{time}}{integer survival time from pulmonary metastasis to death in months} 
#'}
#' @details See Collett (2023)
"pulmonary"

#' @title Recurrence free survival in breast cancer patients
#' @description Recurrence free survival in breast cancer patients
#' @format A data frame with 686 rows and 11 variables:
#' \describe{
#'   \item{\code{id}}{integer patient id}
#'   \item{\code{treat}}{integer hormonal treatment (0=no tamoxifen, 1=tamoxifen)}
#'   \item{\code{age}}{integer patient age (years)}
#'   \item{\code{men}}{integer menopausal status (1=premenopausal, 2=postmenopausal)}
#'   \item{\code{size}}{integer tumour size (mm)}
#'   \item{\code{grade}}{integer tumour grade (1,2,3)}
#'   \item{\code{nodes}}{integer number of positive pymph nodes}
#'   \item{\code{prog}}{integer progesterone receptor status (femtomoles)}
#'   \item{\code{oest}}{integer oestrogen receptor status (femtomoles)}
#'   \item{\code{time}}{integer recurrence-free survival time (days)}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=relapse or death)} 
#'}
#' @details See Collett (2023)
"breast_rfs"

#' @title Recurrence of an ulcer
#' @description A double-blind trial comparing two treatments for ulcers. Data from Belgium.
#' @format A data frame with 43 rows and 6 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{age}}{integer age at the end of the trial in years}
#'   \item{\code{duration}}{integer duration of verified disease (1: <5 years, 2: >=5 years}
#'   \item{\code{treatment}}{integer treatment arm (1=A,2=B) }
#'   \item{\code{time}}{integer time since last visit (months)}
#'   \item{\code{result}}{integer result of the last visit (1=no ulcer detected, 2=ulcer detected)} 
#'}
#' @details See Collett (2023)
"ulcer"

#' @title Recurrence of bladder cancer
#' @description Placebo controlled trial of bladder cancer patients randomised to thiopeta or to placebo
#' @format A data frame with 86 rows and 6 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient number (1-86)}
#'   \item{\code{time}}{integer survival time in months}
#'   \item{\code{status}}{integer status of patient (0=censored, 1=recurrence)}
#'   \item{\code{treat}}{integer treatment group (1=placebo, 2=thiotepa)}
#'   \item{\code{init}}{integer initial number of tumours}
#'   \item{\code{size}}{integer diameter of larger initial tumour in cm} 
#'}
#' @details See Collett (2023)
"bladder"

#' @title Recurrence of mammary tumours in female rats
#' @description This is an animal experiment to compare the use of retinyl acetate (related to vitamin A) across the study (treatment) to treatment with retinyl acetate to 60 days and then no further treatment (control). The female rats all had mammary tumours.
#' @format A data frame with 254 rows and 4 variables:
#' \describe{
#'   \item{\code{rat}}{integer id for each rat}
#'   \item{\code{treatment}}{integer treatment arm indicator (1=treatment, 0=control)}
#'   \item{\code{time}}{double follow-up time (days)}
#'   \item{\code{status}}{integer recurrence indicator (0=no, 1=yes)} 
#'}
#' @details See Collett (2023)
"mammary"

#' @title Survival following aortic valve replacement
#' @description Patients following an aortic valve replacement are measured for left ventricular mass index (LVMI).
#' @format A data frame with 988 rows and 11 variables:
#' \describe{
#'   \item{\code{id}}{integer patient id}
#'   \item{\code{futime}}{double total follow-up time from date of surgery (years)}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=death)}
#'   \item{\code{time}}{double time of LVMI measurement after surgery (years)}
#'   \item{\code{lvmi}}{double standardised LVMI}
#'   \item{\code{age}}{integer age of patient in years}
#'   \item{\code{sex}}{integer sex of patient (0=male, 1=female)}
#'   \item{\code{redo}}{integer previous cardiac surgery (0=no, 1=yes)}
#'   \item{\code{emerg}}{integer operative urgency (0=elective, 1=urgent or emergency)}
#'   \item{\code{dm}}{integer preoperative diabetes mellitus (0=no, 1=yes)}
#'   \item{\code{type}}{integer type of valve (1=human tissue, 2=porcine tissue)} 
#'}
#' @details See Collett (2023)
"valve"

#' @title Survival following kidney transplantation
#' @description Survival following kidney transplantation
#' @format A data frame with 434 rows and 7 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{donor}}{integer donoe id}
#'   \item{\code{time}}{integer survival time in days}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=graft failure or death with a functioning graft)}
#'   \item{\code{age}}{integer patient age (years)}
#'   \item{\code{diabetes}}{integer diabetes status (0=absent, 1=present)}
#'   \item{\code{cit}}{double cold ischaemic time, the time in hours between retrieval of the kidney from the donor and the transplantation} 
#'}
#' @details See Collett (2023)
"tplant"

#' @title Survival of black ducks
#' @description Black ducks, Anas rubripes, were followed the US Fish and Wildlife Service.
#' @format A data frame with 50 rows and 6 variables:
#' \describe{
#'   \item{\code{duck}}{integer duck indicator}
#'   \item{\code{time}}{integer survival time in days}
#'   \item{\code{status}}{integer status of bird (0=alive or missing, 1=dead)}
#'   \item{\code{age}}{integer age group (0=hatch-year bird, 1=bird aged >= 1 year)}
#'   \item{\code{weight}}{integer weight of bird in g}
#'   \item{\code{length}}{integer length of wing in mm} 
#'}
#' @details See Collett (2023)
"ducks"

#' @title Survival of laboratory mice
#' @description Laboratory study of survival for two groups of mice exposed to radiation.
#' @format A data frame with 181 rows and 3 variables:
#' \describe{
#'   \item{\code{environment}}{integer type of environment (1=standard, 2=germ-free)}
#'   \item{\code{causeofdeath}}{integer cause of death (1=thymic lymphoma, 2=reticulum cell sarcoma, 3=other causes)}
#'   \item{\code{time}}{integer survival time (days)} 
#'}
#' @details See Collett (2023). Note that are no censored event times.
"mice"

#' @title Survival of liver transplant recipients
#' @description Survival of liver transplant recipients
#' @format A data frame with 1761 rows and 7 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{age}}{integer patient age in years}
#'   \item{\code{gender}}{integer patient gender (1=male, 2=female)}
#'   \item{\code{disease}}{integer primary disease (1=PBC, 2=PSC, 3=ALD)}
#'   \item{\code{time}}{integer time to event (days)}
#'   \item{\code{status}}{integer cof>0}
#'   \item{\code{cof}}{integer cause of graft failure (0=functioning graft, 1=rejection, 2=thrombosis, 3=recurrent disease, 4=other)} 
#'}
#' @details See Collett (2023)
"liver"

#' @title Survival of multiple myeloma patients
#' @description Patients diagnosed with multiple myeloma who were diagnosed and treated with alkylating agents at West Virginia University Medical Center for ages 50-80 years.
#' @format A data frame with 48 rows and 10 variables:
#' \describe{
#'   \item{\code{patient}}{integer for a patient identifier}
#'   \item{\code{time}}{integer survival time in months}
#'   \item{\code{status}}{integer for status at follow-up (0=Alive, 1=Dead)}
#'   \item{\code{age}}{integer age at diagnosis in years}
#'   \item{\code{sex}}{integer for sex of the patient (1=male, 2=female)}
#'   \item{\code{bun}}{integer level of blood urea nitrogen at diagnosis (unit assumed to be mg/dL based on the normal range for adults reported by \url{https://en.wikipedia.org/wiki/Blood_urea_nitrogen})}
#'   \item{\code{ca}}{integer serum calcium at diagnosis in mg/dL}
#'   \item{\code{hb}}{double for serum hemoglobin level at diagnosis in g/dL (equivalently, grams per 100 mL)}
#'   \item{\code{pcells}}{integer percent of plasma cells in the bone marrow at diagnosis}
#'   \item{\code{protein}}{integer indicator for whether or not the Bence-Jones protein was present in the urine at diagnosis (0=absent, 1=present)} 
#'}
#' @details
#' Krall et al (1975) did not provide the units for all of these measurements. In their analyses, they used some data transformations: log(bun). Collett (2023) converted data from Krall et al (1975): BUN is reported by Krall and colleagues as X1=log(BUN), however the log base and unit is unclear; Krall and colleagues reported for 65 individuals, including those younger than 50 and older than 80.
#'
#' @examples
#' ## To be completed.
#'
#' @references
#' Krall JM, Uthoff VA, Harley JB. A step-up procedure for selecting variables associated with survival. Biometrics. 1975 Mar 1:49-57.
#' \doi{10.2307/2529709}
"myeloma"

#' @title Survival of patients registered for a lung transplant
#' @description Survival of patients registered for a lung transplant
#' @format A data frame with 196 rows and 7 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{time}}{integer time from registration to the earlist of removal from list, last known follow-up date, 30 April 2012, or death (days)}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=dead)}
#'   \item{\code{age}}{integer age in years}
#'   \item{\code{gender}}{integer gender (1=male, 2=female)}
#'   \item{\code{bmi}}{double body mass index}
#'   \item{\code{disease}}{integer disease (1=COPD, 2=fibrosis, 3=suppurative, 4=other)} 
#'}
#' @details See Collett (2023)
"lung"

#' @title Survival of patients with gastric cancer
#' @description Survival of patients with gastric cancer
#' @format A data frame with 90 rows and 4 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{time}}{integer survival time in days}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=dead)}
#'   \item{\code{treat}}{integer treatment arm (0=chemotherapy alone, 1=chemotherapy and radiotherapy)} 
#'}
#' @details See Collett (2023)
"gcancer"

#' @title Survival times of patients with melanoma
#' @description Comparing two immunotherapy treatments for patients with melanoma
#' @format A data frame with 30 rows and 4 variables:
#' \describe{
#'   \item{\code{age}}{integer age group (1=21-44, 2=41-60, 3=61+)}
#'   \item{\code{treatment}}{integer treatment arm (1=BCG, 2=C. parvum)}
#'   \item{\code{time}}{integer survival time (months)}
#'   \item{\code{status}}{integer event indicator (0=censored, 1=dead)} 
#'}
#' @details See Collett (2023)
"melanoma"

#' @title Time to death while waiting for a liver transplant
#' @description Investigate the time on the liver transplantation list.
#' @format A data frame with 281 rows and 7 variables:
#' \describe{
#'   \item{\code{patient}}{integer patient id}
#'   \item{\code{time}}{integer time on the list}
#'   \item{\code{status}}{integer event indicator (0=censored, including having a transplant, 1=died on the list)}
#'   \item{\code{age}}{integer patient age in years}
#'   \item{\code{gender}}{integer patient gender (1=male, 0=female)}
#'   \item{\code{bmi}}{double body mass index (kg/m^2)}
#'   \item{\code{ukeld}}{integer UK endstage liver disease score} 
#'}
#' @details See Collett (2023). A higher UKELD is associated with worse disease severity.
"livertx"

#' @title Time to discontinuation of the use of an IUD
#' @description A very simple dataset showing potential right censoring for time to discontinuation of the use of an IUD.
#' @format A data frame with 18 rows and 2 variables:
#' \describe{
#'   \item{\code{time}}{integer Time in weeks to discontinuation of the use of an IUD}
#'   \item{\code{status}}{integer Indicator for whether the IUD was discontinued: 0=No, 1=Yes} 
#'}
#' @details
#' These data are reported in Table 1.1 (Collett, 2023, page 6). 
"IUD"

#' @title Treatment of hypernephroma
#' @description This study was undertaken at the University of Oklahoma
#'     Health Sciences Center to investigate survival among 36
#'     patients with a kidney tumour (hypernephroma). Standard tangent
#'     included chemotherapy and immunotherapy, with some patients
#'     also having a nephrectomy, or surgical removal of the
#'     kidney. For further details, see Lee and Wang (2013).
#' @format A data frame with 36 rows and 4 variables:
#' \describe{
#'   \item{\code{nephrectomy}}{integer indicator for nephrectomy (0=No; 1=Yes)}
#'   \item{\code{age}}{integer age group (1=<60; 2=60-70; 3=>70}
#'   \item{\code{time}}{integer for the follow-up time in months}
#'   \item{\code{status}}{integer for status at the end of follow-up (1=died; 0=censored)} 
#'}
# ' @details 
#' @references Lee ET, Wang J. Statistical Methods for Survival Data Analysis. New York, NY: John Wiley & Sons; 2013, fourth edition. \url{https://www.wiley.com/en-sg/Statistical+Methods+for+Survival+Data+Analysis%252C+4th+Edition-p-9781118095027}
"kidney"

#' @rdname Datasets
#' @description As an alternative to using the R datasets, the `collett_data`
#'     function allows for reading from the original .dat files that
#'     are stored in the package.
#' @param name Character string with the original filename
#' @examples
#' head(collett_data("A numerical illustration.dat"))
#' ## which is equivalent to: head(illustration)
#' @return A data-frame
#' @export
collett_data <- function(name) {
    new_name = gsub(" ","_",name)
    utils::read.table(system.file(paste0("extdata/",new_name),package="collett"),
                      header=TRUE)
}
