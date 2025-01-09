## Code used to extract and document the datasets
## @source: https://s3-eu-west-1.amazonaws.com/s3-euw1-ap-pe-ws4-cws-documents.ri-prod/9781032252858/Data%20sets%20from%20Modelling%20Survival%20Data%20in%20Medical%20Research%2C%204th%20edition.zip
## The zip file was unzipped into the inst/extdata directory

## For detailed documentation, see R/data.R
## Would this lookup table be helpful for programming?

## meta_data = 
##     list(list(file="A numerical illustration.dat",
##               name="illustration"), # 37*4
##          list(file="Bone marrow transplantation in the treatment of leukaemia.dat",
##               name="leukaemia"), # 23*8
##          list(file="Bone marrow transplantation.dat",
##               name="bone_marrow"), # 37*9
##          list(file="Chemotherapy in ovarian cancer patients.dat",
##               name="ovarian"),  # 26*3
##          list(file="Chronic active hepatitis.dat",
##               name="active_hepatitis"), # 44*3
##          list(file="Chronic granulomatous disease.dat",
##               name="granulomatous"), # 128*12
##          list(file="Clinical trial of tamoxifen in breast cancer patients.dat",
##               name="tamoxifen"), # 641*18
##          list(file="Comparison of two treatments for prostatic cancer.dat",
##               name="prostatic"), # Table 1.4, page 10, Example 1.4, 38*8
##          list(file="Comparisons between kidney transplant centres.dat",
##               name="kidneytx"), # 1439*9
##          list(file="Data from a cirrhosis study (baseline).dat",
##               name="liverbase"), # 12*6
##          list(file="Data from a cirrhosis study (in counting process format).dat",
##               name="liver_counting"), # 54*7
##          list(file="Data from a cirrhosis study (lbr data).dat",
##               name="lbrdata0"), # 42*3
##          list(file="Health evaluation and linkage to primary care.dat",
##               name="HELP"), # 447*7 (name="help" leads to R issues:()
##          list(file="Infection in patients on dialysis.dat",
##               name="dialysis"), # 13*5
##          list(file="Patient outcome following bone marrow transplantation.dat",
##               name="bone_marrow_tx"), # 2204*9
##          list(file="Prognosis for women with breast cancer.dat",
##               name="bcancer"), # Table 1.2, page 7, Example 1.2, 45*3
##          list(file="Pulmonary metastasis.dat",
##               name="pulmonary"), # Example 2.1, page 16, 11*1
##          list(file="Recurrence free survival in breast cancer patients.dat",
##               name="breast_rfs"), # 686*11 
##          list(file="Recurrence of an ulcer.dat",
##               name="ulcer"), # 43*6
##          list(file="Recurrence of bladder cancer.dat",
##               name="bladder"), # 86*6
##          list(file="Recurrence of mammary tumours in female rats.dat",
##               name="mammary"), # 254*4
##          list(file="Survival following aortic valve replacement.dat",
##               name="valve"), # 988*11
##          list(file="Survival following kidney transplantation.dat",
##               name="tplant"), # 434*7
##          list(file="Survival of black ducks.dat",
##               name="ducks"), # 50*6
##          list(file="Survival of laboratory mice.dat",
##               name="mice"), # 181*3
##          list(file="Survival of liver transplant recipients.dat",
##               name="liver"), # 1761*7
##          list(file="Survival of multiple myeloma patients.dat",
##               name="myeloma"), # Table 1.3, page 9, Example 1.3, 48*10
##          list(file="Survival of patients registered for a lung transplant.dat",
##               name="lung"), # 196*7
##          list(file="Survival of patients with gastric cancer.dat",
##               name="gcancer"), # 90*4
##          list(file="Survival times of patients with melanoma .dat",
##               name="melanoma"), # 30*4
##          list(file="Time to death while waiting for a liver transplant.dat",
##               name="livertx"), # 281*7
##          list(file="Time to discontinuation of the use of an IUD.dat",
##               name="IUD"), # Table 1.1, page 6, Example 1.1, 18*2
##          list(file="Treatment of hypernephroma.dat",
##               name="kidney") # 36*4
##          )

## ## Create the R datasets and generate a template for the documentation (#pragma once)
## ## library(collett)
## ## file_names = dir(system.file("extdata",package="collett"), "*.dat")
## library(sinew)
## src_path = "~/src/R/collett"
## file_names = dir(paste0(src_path,"/inst/extdata"), "*.dat")
## ## check that the meta_data file names are complete
## stopifnot(all(sapply(meta_data, "[[", "file") == file_names))
## for (item in meta_data) {
##     assign(item$name, read.table(system.file(paste0("extdata/",item$file),package="collett"),header=TRUE))
##     cat("\n")
##     eval(substitute(sinew::makeOxygen(value), list(value=item$name)))
##     eval(substitute(save(value, file=paste0(src_path, "/data/", item$name, ".rda")),
##                     list(value=item$name)))
## }

## ## Tables for documentation (#pragma once)
## cat("#' \\tabular{ll}{\n")
## cat("#' \\bold{Dataset name} \\tab \\bold{File name} \\cr\n")
## cat("#' -------------------- \\tab ----------------- \\cr\n")
## for (item in meta_data) {
##     cat(sprintf("#' \\link{%s} \\tab \"%s\" \\cr\n", item$name, item$file))
## }
## cat("#' }\n")
##
## cat("#' \\tabular{ll}{\n")
## cat("#' \\bold{Dataset name} \\tab \\bold{File name} \\cr\n")
## cat("#' -------------------- \\tab ----------------- \\cr\n")
## for (item in meta_data[order(sapply(meta_data,"[[","name"))]) {
##     cat(sprintf("#' \\link{%s} \\tab \"%s\" \\cr\n", item$name, item$file))
## }
## cat("#' }\n")
