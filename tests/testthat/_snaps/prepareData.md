# Data preparation with groups (median normalization)

    Code
      prepData[["D"]]
    Output
        state1_1 state1_2 state1_3 state2_1 state2_2 state2_3 state3_1 state3_2
      1       NA       NA       NA       NA       NA       NA       NA       NA
      2       NA       NA       NA       NA       NA       NA       NA       NA
      3       NA       NA       NA       NA       NA       NA       NA       NA
      4       NA       NA       NA       NA       NA       NA       NA       NA
      5       NA       NA       NA       NA       NA       NA       NA       NA
      6       NA       NA       NA       NA       NA       NA       NA       NA
      7       NA       NA       NA       NA       NA       NA       NA       NA
      8       NA       NA       NA       NA       NA       NA       NA       NA
      9       NA       NA       NA       NA       NA       NA       NA       NA
        state3_3
      1       NA
      2       NA
      3       NA
      4       NA
      5       NA
      6       NA
      7       NA
      8       NA
      9       NA

---

    Code
      prepData[["ID"]]
    Output
        peptides proteins
      1     pep1     pro1
      2     pep2     pro2
      3     pep3     pro1
      4     pep4     pro3
      5     pep5     pro3
      6     pep6     pro4
      7     pep7     pro5
      8     pep8     pro2
      9     pep9     pro5

---

    Code
      prepData[["D_long"]]
    Output
      # A tibble: 81 x 3
         name     value group 
         <chr>    <dbl> <fct> 
       1 state1_1    NA state1
       2 state1_2    NA state1
       3 state1_3    NA state1
       4 state2_1    NA state2
       5 state2_2    NA state2
       6 state2_3    NA state2
       7 state3_1    NA state3
       8 state3_2    NA state3
       9 state3_3    NA state3
      10 state1_1    NA state1
      # i 71 more rows

# Data preparation without groups (loess normalization)

    Code
      prepData[["D"]]
    Output
        state1_1 state1_2 state1_3 state2_1 state2_2 state2_3 state3_1 state3_2
      1 26.08083 26.13939 26.41706 26.24908 26.33195       NA 26.15211       NA
      2 22.44117 23.14236 22.52429 22.44392 22.45356       NA 22.44587 22.44373
      3 24.44280 24.42117 24.49493       NA 24.41114       NA 24.41317 24.41317
      4 24.83760 25.01244 24.86740 24.95982 24.96383       NA 25.07829 25.01842
      5 25.67963       NA 25.72126 25.67786 25.67835       NA 25.67835 25.67835
      6 26.24328 26.23675 26.10884 26.29706 26.11750       NA 26.30384 26.24071
      7       NA 23.04771 22.96882 23.25897       NA       NA       NA       NA
      8 26.46535 26.46503 26.43665 26.45803 26.46073       NA 26.46073 26.45301
      9 27.96846 27.96636 27.97211 27.96879 27.96878       NA 27.96878 27.96902
        state3_3
      1 26.33669
      2 21.69832
      3 24.41317
      4       NA
      5 25.67328
      6 26.22603
      7 23.49676
      8 26.43449
      9 27.96967

---

    Code
      prepData[["ID"]]
    Output
        peptides proteins
      1     pep1     pro1
      2     pep2     pro2
      3     pep3     pro1
      4     pep4     pro3
      5     pep5     pro3
      6     pep6     pro4
      7     pep7     pro5
      8     pep8     pro2
      9     pep9     pro5

---

    Code
      prepData[["D_long"]]
    Output
      # A tibble: 81 x 3
         name     value group
         <chr>    <dbl> <lgl>
       1 state1_1  26.1 NA   
       2 state1_2  26.1 NA   
       3 state1_3  26.4 NA   
       4 state2_1  26.2 NA   
       5 state2_2  26.3 NA   
       6 state2_3  NA   NA   
       7 state3_1  26.2 NA   
       8 state3_2  NA   NA   
       9 state3_3  26.3 NA   
      10 state1_1  22.4 NA   
      # i 71 more rows

