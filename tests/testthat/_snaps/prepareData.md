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
      1  25.12073       NA 25.11159 25.33946 25.36341 25.24092 25.33947 25.23668
      2  27.83321       NA 27.81856 27.82100 27.75653 27.85408 27.85990 28.06860
      3  24.65740       NA 24.66017 24.46292 24.63633 26.14245 24.65920 24.71610
      4        NA       NA       NA       NA       NA       NA       NA       NA
      5        NA       NA 24.85677       NA 24.59328       NA       NA       NA
      6  24.58737       NA 24.58466 24.47627 24.26970 24.46738 24.24464 24.71714
      7  25.16835       NA 25.30332 25.20466 25.25371 25.27710       NA 25.21420
      8  25.24499       NA 25.07437 25.03418 25.19411 24.99359 25.04899 25.14118
      9  31.29039       NA 31.28198 31.28994 31.30106 31.27694 31.26804 31.26639
      10 20.76279       NA       NA 20.77217 20.78799 20.32014 20.86961 20.78418
      11 25.71019       NA 25.71717 25.52652 25.54378 25.59879 25.58558 25.60943
      12       NA       NA 22.76986 23.37191       NA 22.95057 22.86325       NA
      13 25.63062       NA 25.65834 25.65957 25.60386 25.52984 25.59256 25.52828
      14 24.05278       NA       NA       NA 22.60754 23.64752       NA       NA
      15 25.76101       NA 25.91484 25.85553 25.95322 25.89927 25.89008 25.87726
         state3_3
      1  25.28443
      2  28.07501
      3  24.66469
      4        NA
      5  22.90831
      6  24.40825
      7  25.14623
      8  25.07151
      9  31.20008
      10 20.86842
      11 25.69762
      12 22.85530
      13 25.39154
      14 23.56779
      15 25.87044

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
       1 state1_1  25.1 NA    NA    
       2 state1_2  NA   NA    NA    
       3 state1_3  25.1 NA    NA    
       4 state2_1  25.3 NA    NA    
       5 state2_2  25.4 NA    NA    
       6 state2_3  25.2 NA    NA    
       7 state3_1  25.3 NA    NA    
       8 state3_2  25.2 NA    NA    
       9 state3_3  25.3 NA    NA    
      10 state1_1  27.8 NA    NA    
      # i 125 more rows

