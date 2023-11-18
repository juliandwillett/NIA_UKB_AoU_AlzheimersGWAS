# For this paper, AoU had egress issues even after filing a request where large files could not be exported without crashing/failing
# I dealt with this by making subsets of each chromosome to be under 100 Mb, which worked better
# You also have to slightly adjust the code to ensure that more than 20 people for each hit, per AoU privacy requirements

for (c in 1:22) {
    # Read in data
    print(glue("Curr chr: {c}"))
    outcome = "AD_any"
    chr = glue("chr{c}")
    data = vroom(glue("aou_step2_rg_{chr}_firthallvariants_{outcome}.regenie"))

    # Ensure exported data meets privacy requirements
    data_min20 = data %>% filter(A1FREQ * N >= 20)
    min(data_min20$A1FREQ * data_min20$N)

    # break into pieces
    data_min20_categories = data_min20 %>% mutate(Category = ifelse(row_number() < 1e6,"A",
                                                               ifelse(row_number() < 2e6,"B",
                                                                      ifelse(row_number()<3e6,"C",
                                                                            ifelse(row_number()<4e6,"D",
                                                                                  ifelse(row_number()<5e6,"E","F"))))))
    subsets = split(data_min20_categories, data_min20_categories$Category)
    
    # write each piece
    for (category in names(subsets)) {
          subset_data <- subsets[[category]]
          filename <- glue("aou_step2_rg_{chr}_allvar_anc_all_{outcome}_min20N_{category}.regenie")  # Define the file name
          vroom_write(subset_data, filename)
    }
}
