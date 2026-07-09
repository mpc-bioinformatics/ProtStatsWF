# ttest returns consistent dimensions and column names

    Code
      dim(ttest_res)
    Output
      [1] 4185   91

---

    Code
      colnames(ttest_res)
    Output
       [1] "Protein"             "Gene"                "Protein.Length"     
       [4] "Organism"            "Description"         "C_1"                
       [7] "C_2"                 "C_3"                 "C_4"                
      [10] "C_5"                 "C_6"                 "C_7"                
      [13] "C_8"                 "C_9"                 "C_10"               
      [16] "C_11"                "C_12"                "C_13"               
      [19] "C_14"                "C_15"                "C_16"               
      [22] "C_17"                "C_18"                "C_19"               
      [25] "HCC_1"               "HCC_2"               "HCC_3"              
      [28] "HCC_4"               "HCC_5"               "HCC_6"              
      [31] "HCC_7"               "HCC_8"               "HCC_9"              
      [34] "HCC_10"              "HCC_11"              "HCC_12"             
      [37] "HCC_13"              "HCC_14"              "HCC_15"             
      [40] "HCC_16"              "HCC_17"              "HCC_18"             
      [43] "HCC_19"              "HCC_1_delog"         "HCC_2_delog"        
      [46] "HCC_3_delog"         "HCC_4_delog"         "HCC_5_delog"        
      [49] "HCC_6_delog"         "HCC_7_delog"         "HCC_8_delog"        
      [52] "HCC_9_delog"         "HCC_10_delog"        "HCC_11_delog"       
      [55] "HCC_12_delog"        "HCC_13_delog"        "HCC_14_delog"       
      [58] "HCC_15_delog"        "HCC_16_delog"        "HCC_17_delog"       
      [61] "HCC_18_delog"        "HCC_19_delog"        "C_1_delog"          
      [64] "C_2_delog"           "C_3_delog"           "C_4_delog"          
      [67] "C_5_delog"           "C_6_delog"           "C_7_delog"          
      [70] "C_8_delog"           "C_9_delog"           "C_10_delog"         
      [73] "C_11_delog"          "C_12_delog"          "C_13_delog"         
      [76] "C_14_delog"          "C_15_delog"          "C_16_delog"         
      [79] "C_17_delog"          "C_18_delog"          "C_19_delog"         
      [82] "mean_of_differences" "test_statistic"      "p"                  
      [85] "p.fdr"               "FC_HCC_divided_by_C" "FC_C_divided_by_HCC"
      [88] "CI_lower"            "CI_upper"            "n_C"                
      [91] "n_HCC"              

# ttest statistical results are consistent

    Code
      ttest_res[1:5, c("p", "p.fdr", fc_col)]
    Output
                         p     p.fdr FC_C_divided_by_HCC
      P00711            NA        NA                  NA
      P00761     0.9426601 0.9648075           1.0017607
      P02666            NA        NA                  NA
      P02769     0.8602800 0.9098981           0.3098479
      A0A075B6I0 0.8862280 0.9297831           0.7732134

