# paired ttest statistical results are consistent

    Code
      ttest_res
    Output
      DataFrame with 4185 rows and 53 columns
                     Protein        Gene Protein.Length     Organism
                 <character> <character>      <integer>  <character>
      P00711          P00711       LALBA            142   Bos taurus
      P00761          P00761          NA            231   Sus scrofa
      P02666          P02666        CSN2            224   Bos taurus
      P02769          P02769         ALB            607   Bos taurus
      A0A075B6I0  A0A075B6I0    IGLV8-61            122 Homo sapiens
      ...                ...         ...            ...          ...
      Q9Y6N7          Q9Y6N7       ROBO1           1651 Homo sapiens
      Q9Y6W5          Q9Y6W5       WASF2            498 Homo sapiens
      Q9Y6X4          Q9Y6X4     FAM169A            670 Homo sapiens
      Q9Y6X9          Q9Y6X9       MORC2           1032 Homo sapiens
      Q9Y6Y8          Q9Y6Y8     SEC23IP           1000 Homo sapiens
                            Description mean_of_differences test_statistic          p
                            <character>           <numeric>      <numeric>  <numeric>
      P00711          Alpha-lactalbumin                  NA             NA         NA
      P00761                    Trypsin          -0.0153754     -0.0729376   0.942660
      P02666                Beta-casein                  NA             NA         NA
      P02769                    Albumin           0.1490828      0.1792886   0.860280
      A0A075B6I0 Immunoglobulin lambd..           0.0999398      0.1483784   0.886228
      ...                           ...                 ...            ...        ...
      Q9Y6N7       Roundabout homolog 1                  NA             NA         NA
      Q9Y6W5     Actin-binding protei..             1.05863        3.89751 0.00115794
      Q9Y6X4     Soluble lamin-associ..                  NA             NA         NA
      Q9Y6X9               ATPase MORC2                  NA             NA         NA
      Q9Y6Y8     SEC23-interacting pr..             0.79051        2.71077 0.01432028
                     p.fdr FC_HCC_divided_by_C FC_C_divided_by_HCC  CI_lower
                 <numeric>           <numeric>           <numeric> <numeric>
      P00711            NA                  NA                  NA        NA
      P00761      0.964807            0.998242            1.001761 -0.458255
      P02666            NA                  NA                  NA        NA
      P02769      0.909898            3.227390            0.309848 -1.634360
      A0A075B6I0  0.929783            1.293304            0.773213 -1.492746
      ...              ...                 ...                 ...       ...
      Q9Y6N7            NA                  NA                  NA        NA
      Q9Y6W5     0.0060150             2.10063            0.476047  0.485568
      Q9Y6X4            NA                  NA                  NA        NA
      Q9Y6X9            NA                  NA                  NA        NA
      Q9Y6Y8     0.0416183             1.43772            0.695544  0.177842
                  CI_upper       n_C     n_HCC       C_1       C_2       C_3
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA   21.0248        NA        NA
      P00761      0.427504        19        19   25.1752   24.4546   24.7485
      P02666            NA        NA        NA        NA        NA        NA
      P02769      1.932525        18        16   23.1098   19.0751   20.3355
      A0A075B6I0  1.692625        11         9   15.2639        NA   18.1863
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA        NA
      Q9Y6W5       1.63169        19        18   17.7218   18.4148   16.7358
      Q9Y6X4            NA        NA        NA        NA        NA        NA
      Q9Y6X9            NA        NA        NA        NA        NA        NA
      Q9Y6Y8       1.40318        19        19   18.8406   19.5081   18.3236
                       C_4       C_5       C_6       C_7       C_8       C_9
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA   19.6445        NA        NA
      P00761       24.0242   24.1876   24.6480   24.7356   23.3475   23.1098
      P02666            NA        NA        NA        NA        NA        NA
      P02769       18.2336   20.8112   21.1691   23.5649   19.4682   22.2155
      A0A075B6I0   17.6219        NA   15.5900   16.2362   17.3138   16.4154
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA        NA
      Q9Y6W5       18.4632   17.8013   18.3259    16.977   18.0377   17.0448
      Q9Y6X4            NA        NA        NA        NA        NA        NA
      Q9Y6X9            NA        NA        NA        NA        NA        NA
      Q9Y6Y8       19.5381   18.7073   16.2949    16.327   19.2177   17.3887
                      C_10      C_11      C_12      C_13      C_14      C_15
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA        NA
      P00761       23.4683   24.0411   24.0662   24.1044   24.0146   24.3000
      P02666            NA        NA        NA        NA        NA        NA
      P02769       20.5234   24.7340   20.6449   17.1282   20.2490   20.7871
      A0A075B6I0   15.4476        NA        NA   17.4574   16.8978        NA
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA        NA
      Q9Y6W5       19.4629   19.0530   18.4927   18.0040   18.9519   17.7989
      Q9Y6X4            NA        NA        NA        NA   17.3557        NA
      Q9Y6X9            NA        NA   18.4978        NA        NA        NA
      Q9Y6Y8       19.8172   19.9588   20.6168   19.3576   20.1798   19.5959
                      C_16      C_17      C_18      C_19     HCC_1     HCC_2
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA   19.3648
      P00761       23.5454   22.7932   22.3884   24.7906   25.1570   24.3500
      P02666            NA        NA        NA        NA        NA        NA
      P02769            NA   19.9835   19.6079   18.3434   18.3731   21.5963
      A0A075B6I0        NA   16.1529        NA        NA        NA        NA
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA   17.6803        NA        NA        NA
      Q9Y6W5       17.0380   17.8139   17.4082   18.6459   19.3255   17.7162
      Q9Y6X4            NA        NA        NA        NA        NA   20.8938
      Q9Y6X9       17.4878   19.2007        NA        NA   17.7774        NA
      Q9Y6Y8       19.3562   19.9587   20.0836   20.5146   19.9889   20.0287
                     HCC_3     HCC_4     HCC_5     HCC_6     HCC_7     HCC_8
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA        NA
      P00761       22.7189   23.0752   23.6687   23.7156   23.6249   23.4768
      P02666            NA        NA        NA        NA        NA        NA
      P02769            NA        NA   21.4039   20.0467        NA   24.0458
      A0A075B6I0   16.1887   15.6806        NA   17.8665        NA   15.8739
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7       18.3014        NA        NA        NA        NA        NA
      Q9Y6W5       19.7710   19.5013   19.2495   19.2590   19.1746   19.8272
      Q9Y6X4       15.4977   22.1396        NA        NA   18.6553   18.9198
      Q9Y6X9       17.1159        NA   18.4473   17.6418        NA        NA
      Q9Y6Y8       20.4184   20.6310   20.1044   19.8790   19.5584   19.5961
                     HCC_9    HCC_10    HCC_11    HCC_12    HCC_13    HCC_14
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA        NA
      P00761       22.3560   24.6541   24.5383   24.1768   23.6484   24.3417
      P02666            NA        NA        NA        NA        NA        NA
      P02769       18.9628   20.4484   20.3110   27.0643   19.7938   23.0814
      A0A075B6I0   17.3053        NA        NA        NA   17.7739   16.5213
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA   18.2904
      Q9Y6W5       19.5974   18.4542        NA   18.9024   19.1817   19.1012
      Q9Y6X4            NA        NA        NA        NA   21.1521        NA
      Q9Y6X9            NA   17.7677   19.5415   19.4161        NA        NA
      Q9Y6Y8       19.4264   20.2950   19.9454   19.3048   19.3752   20.3507
                    HCC_15    HCC_16    HCC_17    HCC_18    HCC_19
                 <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA
      P00761       24.5622   24.6794   24.3745   23.5457   24.9870
      P02666            NA        NA        NA        NA        NA
      P02769       18.1144   21.5357   19.8157   19.3646   17.6642
      A0A075B6I0        NA        NA   19.2248        NA   16.0522
      ...              ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA
      Q9Y6W5       18.4082   17.9691   20.6115   17.7837   18.3604
      Q9Y6X4       17.3789   18.7403   16.5795        NA        NA
      Q9Y6X9       17.9356   17.1239        NA   19.0436   15.8044
      Q9Y6Y8       18.6586   20.1786   19.4617   20.2177   21.1860

# unpaired ttest statistical results are consistent

    Code
      ttest_res
    Output
      DataFrame with 4185 rows and 55 columns
                     Protein        Gene Protein.Length     Organism
                 <character> <character>      <integer>  <character>
      P00711          P00711       LALBA            142   Bos taurus
      P00761          P00761          NA            231   Sus scrofa
      P02666          P02666        CSN2            224   Bos taurus
      P02769          P02769         ALB            607   Bos taurus
      A0A075B6I0  A0A075B6I0    IGLV8-61            122 Homo sapiens
      ...                ...         ...            ...          ...
      Q9Y6N7          Q9Y6N7       ROBO1           1651 Homo sapiens
      Q9Y6W5          Q9Y6W5       WASF2            498 Homo sapiens
      Q9Y6X4          Q9Y6X4     FAM169A            670 Homo sapiens
      Q9Y6X9          Q9Y6X9       MORC2           1032 Homo sapiens
      Q9Y6Y8          Q9Y6Y8     SEC23IP           1000 Homo sapiens
                            Description    mean_C  mean_HCC test_statistic
                            <character> <numeric> <numeric>      <numeric>
      P00711          Alpha-lactalbumin        NA        NA             NA
      P00761                    Trypsin   23.9970   23.9816      -0.063605
      P02666                Beta-casein        NA        NA             NA
      P02769                    Albumin   20.5547   20.7264       0.227117
      A0A075B6I0 Immunoglobulin lambd..   16.5985   16.9430       0.702788
      ...                           ...       ...       ...            ...
      Q9Y6N7       Roundabout homolog 1        NA        NA             NA
      Q9Y6W5     Actin-binding protei..   18.0101   19.0108       3.985338
      Q9Y6X4     Soluble lamin-associ..        NA        NA             NA
      Q9Y6X9               ATPase MORC2   18.3954   17.9650      -0.718736
      Q9Y6Y8     SEC23-interacting pr..   19.1361   19.9266       2.495902
                           p      p.fdr FC_HCC_divided_by_C FC_C_divided_by_HCC
                   <numeric>  <numeric>           <numeric>           <numeric>
      P00711              NA         NA                  NA                  NA
      P00761        0.949637   0.968134                  NA                  NA
      P02666              NA         NA                  NA                  NA
      P02769        0.821946   0.884437                  NA                  NA
      A0A075B6I0    0.492619   0.624555                  NA                  NA
      ...                ...        ...                 ...                 ...
      Q9Y6N7              NA         NA                  NA                  NA
      Q9Y6W5     0.000327716 0.00219056                  NA                  NA
      Q9Y6X4              NA         NA                  NA                  NA
      Q9Y6X9     0.511620877 0.64306131                  NA                  NA
      Q9Y6Y8     0.019493525 0.05381485                  NA                  NA
                  CI_lower  CI_upper       n_C     n_HCC NA_reason_code       C_1
                 <numeric> <numeric> <numeric> <numeric>      <numeric> <numeric>
      P00711            NA        NA        NA        NA              1   21.0248
      P00761     -0.505654  0.474903        19        19             NA   25.1752
      P02666            NA        NA        NA        NA              2        NA
      P02769     -1.375405  1.718828        18        16             NA   23.1098
      A0A075B6I0 -0.697599  1.386719        11         9             NA   15.2639
      ...              ...       ...       ...       ...            ...       ...
      Q9Y6N7            NA        NA        NA        NA              1        NA
      Q9Y6W5      0.490855   1.51053        19        18             NA   17.7218
      Q9Y6X4            NA        NA        NA        NA              1        NA
      Q9Y6X9     -2.085591   1.22475         3        11             NA        NA
      Q9Y6Y8      0.138332   1.44269        19        19             NA   18.8406
                       C_2       C_3       C_4       C_5       C_6       C_7
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA   19.6445
      P00761       24.4546   24.7485   24.0242   24.1876   24.6480   24.7356
      P02666            NA        NA        NA        NA        NA        NA
      P02769       19.0751   20.3355   18.2336   20.8112   21.1691   23.5649
      A0A075B6I0        NA   18.1863   17.6219        NA   15.5900   16.2362
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA        NA
      Q9Y6W5       18.4148   16.7358   18.4632   17.8013   18.3259    16.977
      Q9Y6X4            NA        NA        NA        NA        NA        NA
      Q9Y6X9            NA        NA        NA        NA        NA        NA
      Q9Y6Y8       19.5081   18.3236   19.5381   18.7073   16.2949    16.327
                       C_8       C_9      C_10      C_11      C_12      C_13
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA        NA
      P00761       23.3475   23.1098   23.4683   24.0411   24.0662   24.1044
      P02666            NA        NA        NA        NA        NA        NA
      P02769       19.4682   22.2155   20.5234   24.7340   20.6449   17.1282
      A0A075B6I0   17.3138   16.4154   15.4476        NA        NA   17.4574
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA        NA
      Q9Y6W5       18.0377   17.0448   19.4629   19.0530   18.4927   18.0040
      Q9Y6X4            NA        NA        NA        NA        NA        NA
      Q9Y6X9            NA        NA        NA        NA   18.4978        NA
      Q9Y6Y8       19.2177   17.3887   19.8172   19.9588   20.6168   19.3576
                      C_14      C_15      C_16      C_17      C_18      C_19
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA        NA
      P00761       24.0146   24.3000   23.5454   22.7932   22.3884   24.7906
      P02666            NA        NA        NA        NA        NA        NA
      P02769       20.2490   20.7871        NA   19.9835   19.6079   18.3434
      A0A075B6I0   16.8978        NA        NA   16.1529        NA        NA
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA   17.6803        NA
      Q9Y6W5       18.9519   17.7989   17.0380   17.8139   17.4082   18.6459
      Q9Y6X4       17.3557        NA        NA        NA        NA        NA
      Q9Y6X9            NA        NA   17.4878   19.2007        NA        NA
      Q9Y6Y8       20.1798   19.5959   19.3562   19.9587   20.0836   20.5146
                     HCC_1     HCC_2     HCC_3     HCC_4     HCC_5     HCC_6
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA   19.3648        NA        NA        NA        NA
      P00761       25.1570   24.3500   22.7189   23.0752   23.6687   23.7156
      P02666            NA        NA        NA        NA        NA        NA
      P02769       18.3731   21.5963        NA        NA   21.4039   20.0467
      A0A075B6I0        NA        NA   16.1887   15.6806        NA   17.8665
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA   18.3014        NA        NA        NA
      Q9Y6W5       19.3255   17.7162   19.7710   19.5013   19.2495   19.2590
      Q9Y6X4            NA   20.8938   15.4977   22.1396        NA        NA
      Q9Y6X9       17.7774        NA   17.1159        NA   18.4473   17.6418
      Q9Y6Y8       19.9889   20.0287   20.4184   20.6310   20.1044   19.8790
                     HCC_7     HCC_8     HCC_9    HCC_10    HCC_11    HCC_12
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA        NA
      P00761       23.6249   23.4768   22.3560   24.6541   24.5383   24.1768
      P02666            NA        NA        NA        NA        NA        NA
      P02769            NA   24.0458   18.9628   20.4484   20.3110   27.0643
      A0A075B6I0        NA   15.8739   17.3053        NA        NA        NA
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA        NA        NA        NA        NA        NA
      Q9Y6W5       19.1746   19.8272   19.5974   18.4542        NA   18.9024
      Q9Y6X4       18.6553   18.9198        NA        NA        NA        NA
      Q9Y6X9            NA        NA        NA   17.7677   19.5415   19.4161
      Q9Y6Y8       19.5584   19.5961   19.4264   20.2950   19.9454   19.3048
                    HCC_13    HCC_14    HCC_15    HCC_16    HCC_17    HCC_18
                 <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
      P00711            NA        NA        NA        NA        NA        NA
      P00761       23.6484   24.3417   24.5622   24.6794   24.3745   23.5457
      P02666            NA        NA        NA        NA        NA        NA
      P02769       19.7938   23.0814   18.1144   21.5357   19.8157   19.3646
      A0A075B6I0   17.7739   16.5213        NA        NA   19.2248        NA
      ...              ...       ...       ...       ...       ...       ...
      Q9Y6N7            NA   18.2904        NA        NA        NA        NA
      Q9Y6W5       19.1817   19.1012   18.4082   17.9691   20.6115   17.7837
      Q9Y6X4       21.1521        NA   17.3789   18.7403   16.5795        NA
      Q9Y6X9            NA        NA   17.9356   17.1239        NA   19.0436
      Q9Y6Y8       19.3752   20.3507   18.6586   20.1786   19.4617   20.2177
                    HCC_19
                 <numeric>
      P00711            NA
      P00761       24.9870
      P02666            NA
      P02769       17.6642
      A0A075B6I0   16.0522
      ...              ...
      Q9Y6N7            NA
      Q9Y6W5       18.3604
      Q9Y6X4            NA
      Q9Y6X9       15.8044
      Q9Y6Y8       21.1860

