library(glue)

x = read_csv("GTEx Portal.csv")
s = ""
for (t in unique((x %>% filter(`P-Value` <= 1e-5))$Tissue)) {
    s = paste(s,glue("{t} (p = {min((x %>% filter(Tissue == t))$`P-Value`)})"))
}
