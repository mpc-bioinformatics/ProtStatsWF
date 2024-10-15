# Data preparation with groups (median normalization) for test_file_1

    Code
      prepData[["D"]]
    Output
         state1_1 state1_2 state1_3 state2_1 state2_2 state2_3 state3_1 state3_2
      1        NA       NA       NA       NA       NA       NA       NA       NA
      2        NA       NA       NA       NA       NA       NA       NA       NA
      3        NA       NA       NA       NA       NA       NA       NA       NA
      4        NA       NA       NA       NA       NA       NA       NA       NA
      5        NA       NA       NA       NA       NA       NA       NA       NA
      6        NA       NA       NA       NA       NA       NA       NA       NA
      7        NA       NA       NA       NA       NA       NA       NA       NA
      8        NA       NA       NA       NA       NA       NA       NA       NA
      9        NA       NA       NA       NA       NA       NA       NA       NA
      10       NA       NA       NA       NA       NA       NA       NA       NA
      11       NA       NA       NA       NA       NA       NA       NA       NA
      12       NA       NA       NA       NA       NA       NA       NA       NA
      13       NA       NA       NA       NA       NA       NA       NA       NA
      14       NA       NA       NA       NA       NA       NA       NA       NA
      15       NA       NA       NA       NA       NA       NA       NA       NA
         state3_3
      1        NA
      2        NA
      3        NA
      4        NA
      5        NA
      6        NA
      7        NA
      8        NA
      9        NA
      10       NA
      11       NA
      12       NA
      13       NA
      14       NA
      15       NA

---

    Code
      prepData[["ID"]]
    Output
                      peptide                                protein
      1         AIIEEYLHLNDMK A0A0J9YUS5/E9PVC5/E9PVC6/E9Q9E1/Q6NZJ6
      2  CLAFHDISPQAPTHFLVIPK                          B0R1E3/P70349
      3               EGWEYLK                                 Q8R404
      4             FSLQDPPNK                                 O70475
      5      GPPPTDPYGRPPPYDR            H3BJ30/H3BJW3/H3BKW0/Q6NVF9
      6    KSQVFSTAADGQTQVEIK                                 P38647
      7       LQISHEAAACITALR                      A0A087WPL5/E9QNN1
      8  MLVDDIGDVTITNDGATILK                          F2Z483/P11983
      9              QLIVGVNK            D3YZ68/D3Z3I8/P10126/P62631
      10           RGEDMMHPLK                                 Q9QYJ0
      11          SPLAQMEEERR                   E9Q1G8/E9Q9F5/O55131
      12      TILTLTGVSSLEDVK                                 Q8CHP8
      13      VGEATETALTCLVEK     B1ATS4/B1ATS5/E9Q559/O55143/Q64518
      14              WYLTLAR                                 Q9CRA4
      15         YAALSDQGLDIK                                 Q9ET54

---

    Code
      prepData[["D_long"]]
    Output
      # A tibble: 135 x 4
         name     value group  sample
         <chr>    <dbl> <fct>  <chr> 
       1 state1_1    NA state1 1     
       2 state1_2    NA state1 2     
       3 state1_3    NA state1 3     
       4 state2_1    NA state2 1     
       5 state2_2    NA state2 2     
       6 state2_3    NA state2 3     
       7 state3_1    NA state3 1     
       8 state3_2    NA state3 2     
       9 state3_3    NA state3 3     
      10 state1_1    NA state1 1     
      # i 125 more rows

# Data preparation without groups (loess normalization) for test_file_1

    Code
      prepData[["D"]]
    Output
         state1_1 state1_2 state1_3 state2_1 state2_2 state2_3 state3_1 state3_2
      1  25.17401       NA 25.13042 25.39210 25.30028 25.24947 25.37418 25.25304
      2  27.85402       NA 27.85277 27.85074 27.84997 27.84925 27.85137 27.85151
      3  24.71794       NA 24.69667 24.49887 24.64731 26.18813 24.80531 24.60544
      4        NA       NA       NA       NA       NA       NA       NA       NA
      5        NA       NA 24.88293       NA 25.10851       NA       NA       NA
      6  24.53713       NA 24.52131 24.49102 24.49592 24.45302 24.48873 24.48872
      7  25.24380       NA 25.34360 25.26213 25.19427 25.28867       NA 25.23096
      8  25.38873       NA 25.16839 25.09942 25.14052 25.01702 25.12461 25.15537
      9  31.25218       NA 31.25251 31.25293 31.25321 31.25333 31.25287 31.25272
      10 20.79093       NA       NA 20.78549 20.79299 20.79050 20.78885 20.79058
      11 25.68150       NA 25.61966 25.53219 25.44982 25.56617 25.54611 25.57782
      12       NA       NA 22.98156 22.99282       NA 22.96825 22.98763       NA
      13 25.60970       NA 25.57355 25.67214 25.51303 25.50290 25.56196 25.50430
      14 23.58371       NA       NA       NA 23.58957 23.63101       NA       NA
      15 25.69314       NA 25.75582 25.82146 25.84741 25.83440 25.81049 25.81518
         state3_3
      1  25.30463
      2  27.85021
      3  24.77626
      4        NA
      5  22.90047
      6  24.53020
      7  25.17657
      8  25.13651
      9  31.25243
      10 20.79478
      11 25.64942
      12 22.96895
      13 25.35058
      14 23.62305
      15 25.79201

---

    Code
      prepData[["ID"]]
    Output
                      peptide                                protein
      1         AIIEEYLHLNDMK A0A0J9YUS5/E9PVC5/E9PVC6/E9Q9E1/Q6NZJ6
      2  CLAFHDISPQAPTHFLVIPK                          B0R1E3/P70349
      3               EGWEYLK                                 Q8R404
      4             FSLQDPPNK                                 O70475
      5      GPPPTDPYGRPPPYDR            H3BJ30/H3BJW3/H3BKW0/Q6NVF9
      6    KSQVFSTAADGQTQVEIK                                 P38647
      7       LQISHEAAACITALR                      A0A087WPL5/E9QNN1
      8  MLVDDIGDVTITNDGATILK                          F2Z483/P11983
      9              QLIVGVNK            D3YZ68/D3Z3I8/P10126/P62631
      10           RGEDMMHPLK                                 Q9QYJ0
      11          SPLAQMEEERR                   E9Q1G8/E9Q9F5/O55131
      12      TILTLTGVSSLEDVK                                 Q8CHP8
      13      VGEATETALTCLVEK     B1ATS4/B1ATS5/E9Q559/O55143/Q64518
      14              WYLTLAR                                 Q9CRA4
      15         YAALSDQGLDIK                                 Q9ET54

---

    Code
      prepData[["D_long"]]
    Output
      # A tibble: 135 x 4
         name     value group sample
         <chr>    <dbl> <lgl> <lgl> 
       1 state1_1  25.2 NA    NA    
       2 state1_2  NA   NA    NA    
       3 state1_3  25.1 NA    NA    
       4 state2_1  25.4 NA    NA    
       5 state2_2  25.3 NA    NA    
       6 state2_3  25.2 NA    NA    
       7 state3_1  25.4 NA    NA    
       8 state3_2  25.3 NA    NA    
       9 state3_3  25.3 NA    NA    
      10 state1_1  27.9 NA    NA    
      # i 125 more rows

