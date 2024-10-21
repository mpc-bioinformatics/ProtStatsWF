# Calculate ttest 

    Code
      data[["D"]]
    Output
           state1_1 state1_2   state1_3   state2_1   state2_2   state2_3
      1    34910128       NA   33763808   41557998   43313220   39106931
      2   228724928       NA  230557560  228792416  217238016  233899104
      3    25142570       NA   24140170   22372790   26638020   74038624
      4          NA       NA         NA         NA         NA         NA
      5          NA       NA   26402860         NA   26541960         NA
      6    23729240       NA   22417890   22312750   20940040   23442920
      7    36068288       NA   38500928   37816312   40211580   40147352
      8    37990528       NA   32654220   33514620   38775700   33104830
      9  2590670618       NA 2707158182 2511590350 2452957074 2497670980
      10    1495972       NA         NA    1488327    2110030    1490940
      11   52227550       NA   52181250   47326970   48179732   49541580
      12         NA       NA    5683145    9846832         NA    8586830
      13   49514448       NA   50004128   51948072   50367920   47302128
      14   15940000       NA         NA         NA    6838328   13688880
      15   53367728       NA   60332520   58853440   62984728   60534832

---

    Code
      data[["ID"]]
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
      pData
    Output
                      peptide                                protein   state1_1
      1         AIIEEYLHLNDMK A0A0J9YUS5/E9PVC5/E9PVC6/E9Q9E1/Q6NZJ6   34910128
      2  CLAFHDISPQAPTHFLVIPK                          B0R1E3/P70349  228724928
      3               EGWEYLK                                 Q8R404   25142570
      4             FSLQDPPNK                                 O70475         NA
      5      GPPPTDPYGRPPPYDR            H3BJ30/H3BJW3/H3BKW0/Q6NVF9         NA
      6    KSQVFSTAADGQTQVEIK                                 P38647   23729240
      7       LQISHEAAACITALR                      A0A087WPL5/E9QNN1   36068288
      8  MLVDDIGDVTITNDGATILK                          F2Z483/P11983   37990528
      9              QLIVGVNK            D3YZ68/D3Z3I8/P10126/P62631 2590670618
      10           RGEDMMHPLK                                 Q9QYJ0    1495972
      11          SPLAQMEEERR                   E9Q1G8/E9Q9F5/O55131   52227550
      12      TILTLTGVSSLEDVK                                 Q8CHP8         NA
      13      VGEATETALTCLVEK     B1ATS4/B1ATS5/E9Q559/O55143/Q64518   49514448
      14              WYLTLAR                                 Q9CRA4   15940000
      15         YAALSDQGLDIK                                 Q9ET54   53367728
         state1_2   state1_3   state2_1   state2_2   state2_3 state1_1_log
      1        NA   33763808   41557998   43313220   39106931     25.05714
      2        NA  230557560  228792416  217238016  233899104     27.76904
      3        NA   24140170   22372790   26638020   74038624     24.58363
      4        NA         NA         NA         NA         NA           NA
      5        NA   26402860         NA   26541960         NA           NA
      6        NA   22417890   22312750   20940040   23442920     24.50016
      7        NA   38500928   37816312   40211580   40147352     25.10423
      8        NA   32654220   33514620   38775700   33104830     25.17914
      9        NA 2707158182 2511590350 2452957074 2497670980     31.27068
      10       NA         NA    1488327    2110030    1490940     20.51265
      11       NA   52181250   47326970   48179732   49541580     25.63831
      12       NA    5683145    9846832         NA    8586830           NA
      13       NA   50004128   51948072   50367920   47302128     25.56135
      14       NA         NA         NA    6838328   13688880     23.92615
      15       NA   60332520   58853440   62984728   60534832     25.66946
         state1_2_log state1_3_log state2_1_log state2_2_log state2_3_log mean_state1
      1            NA     25.00897     25.30862     25.36830     25.22092          NA
      2            NA     27.78055     27.76946     27.69470     27.80131          NA
      3            NA     24.52493     24.41524     24.66698     26.14177          NA
      4            NA           NA           NA           NA           NA          NA
      5            NA     24.65419           NA     24.66177           NA          NA
      6            NA     24.41815     24.41136     24.31976     24.48265          NA
      7            NA     25.19839     25.17251     25.26111     25.25880          NA
      8            NA     24.96077     24.99829     25.20865     24.98054          NA
      9            NA     31.33413     31.22595     31.19187     31.21794          NA
      10           NA           NA     20.50526     21.00883     20.50779          NA
      11           NA     25.63703     25.49616     25.52192     25.56214          NA
      12           NA     22.43826     23.23123           NA     23.03369          NA
      13           NA     25.57554     25.63057     25.58600     25.49540          NA
      14           NA           NA           NA     22.70521     23.70650          NA
      15           NA     25.84643     25.81062     25.90850     25.85126          NA
         mean_state2 test_statistic  p p.fdr FC_state2_divided_by_state1
      1           NA             NA NA    NA                          NA
      2           NA             NA NA    NA                          NA
      3           NA             NA NA    NA                          NA
      4           NA             NA NA    NA                          NA
      5           NA             NA NA    NA                          NA
      6           NA             NA NA    NA                          NA
      7           NA             NA NA    NA                          NA
      8           NA             NA NA    NA                          NA
      9           NA             NA NA    NA                          NA
      10          NA             NA NA    NA                          NA
      11          NA             NA NA    NA                          NA
      12          NA             NA NA    NA                          NA
      13          NA             NA NA    NA                          NA
      14          NA             NA NA    NA                          NA
      15          NA             NA NA    NA                          NA
         FC_state1_divided_by_state2 CI_lower CI_upper n_state1 n_state2
      1                           NA       NA       NA       NA       NA
      2                           NA       NA       NA       NA       NA
      3                           NA       NA       NA       NA       NA
      4                           NA       NA       NA       NA       NA
      5                           NA       NA       NA       NA       NA
      6                           NA       NA       NA       NA       NA
      7                           NA       NA       NA       NA       NA
      8                           NA       NA       NA       NA       NA
      9                           NA       NA       NA       NA       NA
      10                          NA       NA       NA       NA       NA
      11                          NA       NA       NA       NA       NA
      12                          NA       NA       NA       NA       NA
      13                          NA       NA       NA       NA       NA
      14                          NA       NA       NA       NA       NA
      15                          NA       NA       NA       NA       NA
         NA_reason_code
      1               1
      2               1
      3               1
      4               2
      5               1
      6               1
      7               1
      8               1
      9               1
      10              1
      11              1
      12              1
      13              1
      14              1
      15              1

