# Data preparation with groups (median normalization) for test_file_1

    Code
      prepData$SE
    Output
      # A SummarizedExperiment-tibble abstraction: 144 x 9
      # Features=16 | Samples=9 | Assays=intensity_norm, intensity
         .feature  .sample intensity_norm intensity SampleName group replicate peptide
         <chr>     <chr>            <dbl>     <dbl> <chr>      <chr>     <dbl> <chr>  
       1 AALNALQP~ state1~           26.7      26.3 state1_1   stat~         1 AALNAL~
       2 ESSSHHPG~ state1~           32.0      31.5 state1_1   stat~         1 ESSSHH~
       3 FEAHPNDL~ state1~           NA        NA   state1_1   stat~         1 FEAHPN~
       4 LPNSVLGR  state1~           25.2      24.8 state1_1   stat~         1 LPNSVL~
       5 LQQDIEAVK state1~           29.3      28.8 state1_1   stat~         1 LQQDIE~
       6 MELERPGG~ state1~           33.3      32.8 state1_1   stat~         1 MELERP~
       7 MQEVVANL~ state1~           22.5      22.1 state1_1   stat~         1 MQEVVA~
       8 NNKDSHSL~ state1~           29.5      29.0 state1_1   stat~         1 NNKDSH~
       9 NVVIAADG~ state1~           25.1      24.7 state1_1   stat~         1 NVVIAA~
      10 QGQDGLLS~ state1~           25.3      24.9 state1_1   stat~         1 QGQDGL~
      # i 134 more rows
      # i 1 more variable: protein <chr>

---

    Code
      SummarizedExperiment::colData(prepData$SE)
    Output
      DataFrame with 9 rows and 3 columns
                SampleName       group replicate
               <character> <character> <numeric>
      state1_1    state1_1      state1         1
      state1_2    state1_2      state1         2
      state1_3    state1_3      state1         3
      state2_1    state2_1      state2         1
      state2_2    state2_2      state2         2
      state2_3    state2_3      state2         3
      state3_1    state3_1      state3         1
      state3_2    state3_2      state3         2
      state3_3    state3_3      state3         3

---

    Code
      SummarizedExperiment::rowData(prepData$SE)
    Output
      DataFrame with 16 rows and 2 columns
                                        peptide                protein
                                    <character>            <character>
      AALNALQPPEFR                 AALNALQPPEFR                 Q78T54
      ESSSHHPGIAEFPSR           ESSSHHPGIAEFPSR                 P02671
      FEAHPNDLYVEGLPENIPFR FEAHPNDLYVEGLPENIPFR   G3UYD0/G3UYJ6/Q3UHU8
      LPNSVLGR                         LPNSVLGR                 Q8BH64
      LQQDIEAVK                       LQQDIEAVK                 Q99M74
      ...                                   ...                    ...
      RLDETPDGRK                     RLDETPDGRK                 Q61586
      VEADIAGHGQEVLIR           VEADIAGHGQEVLIR contaminant_MYG_HORS..
      VGVTVAQTTMEPHLLEACVR VGVTVAQTTMEPHLLEACVR A0A0A6YY08/D3Z158/Q8..
      VVAGVANALAHK                 VVAGVANALAHK contaminant_HBB_HUMA..
      YTPLYPFR                         YTPLYPFR                 Q9DCU6

---

    Code
      SummarizedExperiment::assays(prepData$SE)$intensity
    Output
                           state1_1 state1_2 state1_3 state2_1 state2_2 state2_3
      AALNALQPPEFR         26.29200 26.27750 26.02607 26.27129 26.18768 26.02014
      ESSSHHPGIAEFPSR      31.51382 31.71444 31.71987 30.07252 30.24838 30.26180
      FEAHPNDLYVEGLPENIPFR       NA       NA 24.14676 23.75437       NA 23.67309
      LPNSVLGR             24.84863 24.75289 24.88012 24.79729 24.84595 24.70932
      LQQDIEAVK            28.84429 28.65192 28.78472 28.80756 28.67331 28.72744
      MELERPGGNEITR        32.76103 32.80139 32.83883 31.37907 31.35600 31.22309
      MQEVVANLQYDDGSGMK    22.12973 22.27397       NA       NA 22.27230 22.65332
      NNKDSHSLTTNIMEILR    29.02853 29.23454 29.12969 26.35537 28.04455 27.90173
      NVVIAADGVLK          24.74101       NA       NA 24.70867 24.62406 24.78887
      QGQDGLLSVK           24.94521 25.03858 24.63839 25.01088 25.01545 24.93329
      QHIEKAK              24.16319 24.86101 24.62995 22.39856 23.12866 23.43029
      RLDETPDGRK           23.15218 23.79672 24.02570 21.33880 22.11787 22.24675
      VEADIAGHGQEVLIR      28.16441 28.22830 28.17813 25.56709 25.59122 25.59452
      VGVTVAQTTMEPHLLEACVR 24.02252 24.26503       NA 24.20852 24.28883 24.23009
      VVAGVANALAHK         27.48845 27.26185 27.23927 31.03652 31.14066 31.08577
      YTPLYPFR             22.65456 22.96563 22.73358 24.23959 24.28901 23.58755
                           state3_1 state3_2 state3_3
      AALNALQPPEFR         26.32325 26.12448 25.98013
      ESSSHHPGIAEFPSR      28.00825 28.03848 28.19796
      FEAHPNDLYVEGLPENIPFR 22.38642 24.32594 24.40562
      LPNSVLGR             24.93245 24.68149 24.94754
      LQQDIEAVK            28.88355 28.87534 28.72477
      MELERPGGNEITR        29.12220 29.17937 28.97202
      MQEVVANLQYDDGSGMK    22.17315       NA 22.13595
      NNKDSHSLTTNIMEILR    24.33515 25.54931 25.21031
      NVVIAADGVLK                NA       NA       NA
      QGQDGLLSVK           25.16654 24.94118 25.18352
      QHIEKAK                    NA       NA       NA
      RLDETPDGRK                 NA       NA       NA
      VEADIAGHGQEVLIR      29.73262 29.81801 29.76337
      VGVTVAQTTMEPHLLEACVR 24.17893       NA 24.15765
      VVAGVANALAHK         32.23194 32.20771 32.31219
      YTPLYPFR             24.38244 23.82405 23.53935

---

    Code
      SummarizedExperiment::assays(prepData$SE)$intensity_norm
    Output
                           state1_1 state1_2 state1_3 state2_1 state2_2 state2_3
      AALNALQPPEFR         26.69899 25.94292 25.33135 26.60791 26.51839 26.51234
      ESSSHHPGIAEFPSR      32.00165 31.31065 30.87318 30.45785 30.63038 30.83423
      FEAHPNDLYVEGLPENIPFR       NA       NA 23.50221 24.05874       NA 24.12089
      LPNSVLGR             25.23329 24.43773 24.21600 25.11502 25.15972 25.17672
      LQQDIEAVK            29.29079 28.28711 28.01637 29.17668 29.03541 29.27085
      MELERPGGNEITR        33.26817 32.38376 31.96226 31.78114 31.75198 31.81370
      MQEVVANLQYDDGSGMK    22.47230 21.99037       NA       NA 22.55357 23.08183
      NNKDSHSLTTNIMEILR    29.47789 28.86231 28.35214 26.69306 28.39871 28.42952
      NVVIAADGVLK          25.12400       NA       NA 25.02527 24.93503 25.25777
      QGQDGLLSVK           25.33135 24.71978 23.98072 25.33135 25.33135 25.40493
      QHIEKAK              24.53723 24.54448 23.97250 22.68556 23.42074 23.87350
      RLDETPDGRK           23.51057 23.49373 23.38438 21.61222 22.39719 22.66757
      VEADIAGHGQEVLIR      28.60039 27.86889 27.42597 25.89468 25.91440 26.07866
      VGVTVAQTTMEPHLLEACVR 24.39439 23.95608       NA 24.51871 24.59557 24.68842
      VVAGVANALAHK         27.91397 26.91475 26.51218 31.43420 31.53392 31.67379
      YTPLYPFR             23.00525 22.67322 22.12675 24.55018 24.59575 24.03374
                           state3_1 state3_2 state3_3
      AALNALQPPEFR         26.49564 25.33135 26.10488
      ESSSHHPGIAEFPSR      28.19168 27.18724 28.33335
      FEAHPNDLYVEGLPENIPFR 22.53303 23.58741 24.52281
      LPNSVLGR             25.09574 23.93217 25.06732
      LQQDIEAVK            29.07271 27.99870 28.86269
      MELERPGGNEITR        29.31292 28.29350 29.11113
      MQEVVANLQYDDGSGMK    22.31836       NA 22.24224
      NNKDSHSLTTNIMEILR    24.49452 24.77364 25.33135
      NVVIAADGVLK                NA       NA       NA
      QGQDGLLSVK           25.33135 24.18398 25.30444
      QHIEKAK                    NA       NA       NA
      RLDETPDGRK                 NA       NA       NA
      VEADIAGHGQEVLIR      29.92734 28.91275 29.90628
      VGVTVAQTTMEPHLLEACVR 24.33728       NA 24.27364
      VVAGVANALAHK         32.44303 31.22990 32.46734
      YTPLYPFR             24.54213 23.10076 23.65237

---

    Code
      print(prepData$D_long, n = 200)
    Output
      # A tibble: 144 x 8
          .feature .sample intensity_norm intensity SampleName group replicate protein
          <chr>    <fct>            <dbl>     <dbl> <chr>      <chr>     <dbl> <chr>  
        1 AALNALQ~ state1~           26.7      26.3 state1_1   stat~         1 Q78T54 
        2 ESSSHHP~ state1~           32.0      31.5 state1_1   stat~         1 P02671 
        3 FEAHPND~ state1~           NA        NA   state1_1   stat~         1 G3UYD0~
        4 LPNSVLGR state1~           25.2      24.8 state1_1   stat~         1 Q8BH64 
        5 LQQDIEA~ state1~           29.3      28.8 state1_1   stat~         1 Q99M74 
        6 MELERPG~ state1~           33.3      32.8 state1_1   stat~         1 P02671 
        7 MQEVVAN~ state1~           22.5      22.1 state1_1   stat~         1 Q3U1J4 
        8 NNKDSHS~ state1~           29.5      29.0 state1_1   stat~         1 P02671 
        9 NVVIAAD~ state1~           25.1      24.7 state1_1   stat~         1 Q9ESL4 
       10 QGQDGLL~ state1~           25.3      24.9 state1_1   stat~         1 P57759 
       11 QHIEKAK  state1~           24.5      24.2 state1_1   stat~         1 A0A087~
       12 RLDETPD~ state1~           23.5      23.2 state1_1   stat~         1 Q61586 
       13 VEADIAG~ state1~           28.6      28.2 state1_1   stat~         1 contam~
       14 VGVTVAQ~ state1~           24.4      24.0 state1_1   stat~         1 A0A0A6~
       15 VVAGVAN~ state1~           27.9      27.5 state1_1   stat~         1 contam~
       16 YTPLYPFR state1~           23.0      22.7 state1_1   stat~         1 Q9DCU6 
       17 AALNALQ~ state1~           25.9      26.3 state1_2   stat~         2 Q78T54 
       18 ESSSHHP~ state1~           31.3      31.7 state1_2   stat~         2 P02671 
       19 FEAHPND~ state1~           NA        NA   state1_2   stat~         2 G3UYD0~
       20 LPNSVLGR state1~           24.4      24.8 state1_2   stat~         2 Q8BH64 
       21 LQQDIEA~ state1~           28.3      28.7 state1_2   stat~         2 Q99M74 
       22 MELERPG~ state1~           32.4      32.8 state1_2   stat~         2 P02671 
       23 MQEVVAN~ state1~           22.0      22.3 state1_2   stat~         2 Q3U1J4 
       24 NNKDSHS~ state1~           28.9      29.2 state1_2   stat~         2 P02671 
       25 NVVIAAD~ state1~           NA        NA   state1_2   stat~         2 Q9ESL4 
       26 QGQDGLL~ state1~           24.7      25.0 state1_2   stat~         2 P57759 
       27 QHIEKAK  state1~           24.5      24.9 state1_2   stat~         2 A0A087~
       28 RLDETPD~ state1~           23.5      23.8 state1_2   stat~         2 Q61586 
       29 VEADIAG~ state1~           27.9      28.2 state1_2   stat~         2 contam~
       30 VGVTVAQ~ state1~           24.0      24.3 state1_2   stat~         2 A0A0A6~
       31 VVAGVAN~ state1~           26.9      27.3 state1_2   stat~         2 contam~
       32 YTPLYPFR state1~           22.7      23.0 state1_2   stat~         2 Q9DCU6 
       33 AALNALQ~ state1~           25.3      26.0 state1_3   stat~         3 Q78T54 
       34 ESSSHHP~ state1~           30.9      31.7 state1_3   stat~         3 P02671 
       35 FEAHPND~ state1~           23.5      24.1 state1_3   stat~         3 G3UYD0~
       36 LPNSVLGR state1~           24.2      24.9 state1_3   stat~         3 Q8BH64 
       37 LQQDIEA~ state1~           28.0      28.8 state1_3   stat~         3 Q99M74 
       38 MELERPG~ state1~           32.0      32.8 state1_3   stat~         3 P02671 
       39 MQEVVAN~ state1~           NA        NA   state1_3   stat~         3 Q3U1J4 
       40 NNKDSHS~ state1~           28.4      29.1 state1_3   stat~         3 P02671 
       41 NVVIAAD~ state1~           NA        NA   state1_3   stat~         3 Q9ESL4 
       42 QGQDGLL~ state1~           24.0      24.6 state1_3   stat~         3 P57759 
       43 QHIEKAK  state1~           24.0      24.6 state1_3   stat~         3 A0A087~
       44 RLDETPD~ state1~           23.4      24.0 state1_3   stat~         3 Q61586 
       45 VEADIAG~ state1~           27.4      28.2 state1_3   stat~         3 contam~
       46 VGVTVAQ~ state1~           NA        NA   state1_3   stat~         3 A0A0A6~
       47 VVAGVAN~ state1~           26.5      27.2 state1_3   stat~         3 contam~
       48 YTPLYPFR state1~           22.1      22.7 state1_3   stat~         3 Q9DCU6 
       49 AALNALQ~ state2~           26.6      26.3 state2_1   stat~         1 Q78T54 
       50 ESSSHHP~ state2~           30.5      30.1 state2_1   stat~         1 P02671 
       51 FEAHPND~ state2~           24.1      23.8 state2_1   stat~         1 G3UYD0~
       52 LPNSVLGR state2~           25.1      24.8 state2_1   stat~         1 Q8BH64 
       53 LQQDIEA~ state2~           29.2      28.8 state2_1   stat~         1 Q99M74 
       54 MELERPG~ state2~           31.8      31.4 state2_1   stat~         1 P02671 
       55 MQEVVAN~ state2~           NA        NA   state2_1   stat~         1 Q3U1J4 
       56 NNKDSHS~ state2~           26.7      26.4 state2_1   stat~         1 P02671 
       57 NVVIAAD~ state2~           25.0      24.7 state2_1   stat~         1 Q9ESL4 
       58 QGQDGLL~ state2~           25.3      25.0 state2_1   stat~         1 P57759 
       59 QHIEKAK  state2~           22.7      22.4 state2_1   stat~         1 A0A087~
       60 RLDETPD~ state2~           21.6      21.3 state2_1   stat~         1 Q61586 
       61 VEADIAG~ state2~           25.9      25.6 state2_1   stat~         1 contam~
       62 VGVTVAQ~ state2~           24.5      24.2 state2_1   stat~         1 A0A0A6~
       63 VVAGVAN~ state2~           31.4      31.0 state2_1   stat~         1 contam~
       64 YTPLYPFR state2~           24.6      24.2 state2_1   stat~         1 Q9DCU6 
       65 AALNALQ~ state2~           26.5      26.2 state2_2   stat~         2 Q78T54 
       66 ESSSHHP~ state2~           30.6      30.2 state2_2   stat~         2 P02671 
       67 FEAHPND~ state2~           NA        NA   state2_2   stat~         2 G3UYD0~
       68 LPNSVLGR state2~           25.2      24.8 state2_2   stat~         2 Q8BH64 
       69 LQQDIEA~ state2~           29.0      28.7 state2_2   stat~         2 Q99M74 
       70 MELERPG~ state2~           31.8      31.4 state2_2   stat~         2 P02671 
       71 MQEVVAN~ state2~           22.6      22.3 state2_2   stat~         2 Q3U1J4 
       72 NNKDSHS~ state2~           28.4      28.0 state2_2   stat~         2 P02671 
       73 NVVIAAD~ state2~           24.9      24.6 state2_2   stat~         2 Q9ESL4 
       74 QGQDGLL~ state2~           25.3      25.0 state2_2   stat~         2 P57759 
       75 QHIEKAK  state2~           23.4      23.1 state2_2   stat~         2 A0A087~
       76 RLDETPD~ state2~           22.4      22.1 state2_2   stat~         2 Q61586 
       77 VEADIAG~ state2~           25.9      25.6 state2_2   stat~         2 contam~
       78 VGVTVAQ~ state2~           24.6      24.3 state2_2   stat~         2 A0A0A6~
       79 VVAGVAN~ state2~           31.5      31.1 state2_2   stat~         2 contam~
       80 YTPLYPFR state2~           24.6      24.3 state2_2   stat~         2 Q9DCU6 
       81 AALNALQ~ state2~           26.5      26.0 state2_3   stat~         3 Q78T54 
       82 ESSSHHP~ state2~           30.8      30.3 state2_3   stat~         3 P02671 
       83 FEAHPND~ state2~           24.1      23.7 state2_3   stat~         3 G3UYD0~
       84 LPNSVLGR state2~           25.2      24.7 state2_3   stat~         3 Q8BH64 
       85 LQQDIEA~ state2~           29.3      28.7 state2_3   stat~         3 Q99M74 
       86 MELERPG~ state2~           31.8      31.2 state2_3   stat~         3 P02671 
       87 MQEVVAN~ state2~           23.1      22.7 state2_3   stat~         3 Q3U1J4 
       88 NNKDSHS~ state2~           28.4      27.9 state2_3   stat~         3 P02671 
       89 NVVIAAD~ state2~           25.3      24.8 state2_3   stat~         3 Q9ESL4 
       90 QGQDGLL~ state2~           25.4      24.9 state2_3   stat~         3 P57759 
       91 QHIEKAK  state2~           23.9      23.4 state2_3   stat~         3 A0A087~
       92 RLDETPD~ state2~           22.7      22.2 state2_3   stat~         3 Q61586 
       93 VEADIAG~ state2~           26.1      25.6 state2_3   stat~         3 contam~
       94 VGVTVAQ~ state2~           24.7      24.2 state2_3   stat~         3 A0A0A6~
       95 VVAGVAN~ state2~           31.7      31.1 state2_3   stat~         3 contam~
       96 YTPLYPFR state2~           24.0      23.6 state2_3   stat~         3 Q9DCU6 
       97 AALNALQ~ state3~           26.5      26.3 state3_1   stat~         1 Q78T54 
       98 ESSSHHP~ state3~           28.2      28.0 state3_1   stat~         1 P02671 
       99 FEAHPND~ state3~           22.5      22.4 state3_1   stat~         1 G3UYD0~
      100 LPNSVLGR state3~           25.1      24.9 state3_1   stat~         1 Q8BH64 
      101 LQQDIEA~ state3~           29.1      28.9 state3_1   stat~         1 Q99M74 
      102 MELERPG~ state3~           29.3      29.1 state3_1   stat~         1 P02671 
      103 MQEVVAN~ state3~           22.3      22.2 state3_1   stat~         1 Q3U1J4 
      104 NNKDSHS~ state3~           24.5      24.3 state3_1   stat~         1 P02671 
      105 NVVIAAD~ state3~           NA        NA   state3_1   stat~         1 Q9ESL4 
      106 QGQDGLL~ state3~           25.3      25.2 state3_1   stat~         1 P57759 
      107 QHIEKAK  state3~           NA        NA   state3_1   stat~         1 A0A087~
      108 RLDETPD~ state3~           NA        NA   state3_1   stat~         1 Q61586 
      109 VEADIAG~ state3~           29.9      29.7 state3_1   stat~         1 contam~
      110 VGVTVAQ~ state3~           24.3      24.2 state3_1   stat~         1 A0A0A6~
      111 VVAGVAN~ state3~           32.4      32.2 state3_1   stat~         1 contam~
      112 YTPLYPFR state3~           24.5      24.4 state3_1   stat~         1 Q9DCU6 
      113 AALNALQ~ state3~           25.3      26.1 state3_2   stat~         2 Q78T54 
      114 ESSSHHP~ state3~           27.2      28.0 state3_2   stat~         2 P02671 
      115 FEAHPND~ state3~           23.6      24.3 state3_2   stat~         2 G3UYD0~
      116 LPNSVLGR state3~           23.9      24.7 state3_2   stat~         2 Q8BH64 
      117 LQQDIEA~ state3~           28.0      28.9 state3_2   stat~         2 Q99M74 
      118 MELERPG~ state3~           28.3      29.2 state3_2   stat~         2 P02671 
      119 MQEVVAN~ state3~           NA        NA   state3_2   stat~         2 Q3U1J4 
      120 NNKDSHS~ state3~           24.8      25.5 state3_2   stat~         2 P02671 
      121 NVVIAAD~ state3~           NA        NA   state3_2   stat~         2 Q9ESL4 
      122 QGQDGLL~ state3~           24.2      24.9 state3_2   stat~         2 P57759 
      123 QHIEKAK  state3~           NA        NA   state3_2   stat~         2 A0A087~
      124 RLDETPD~ state3~           NA        NA   state3_2   stat~         2 Q61586 
      125 VEADIAG~ state3~           28.9      29.8 state3_2   stat~         2 contam~
      126 VGVTVAQ~ state3~           NA        NA   state3_2   stat~         2 A0A0A6~
      127 VVAGVAN~ state3~           31.2      32.2 state3_2   stat~         2 contam~
      128 YTPLYPFR state3~           23.1      23.8 state3_2   stat~         2 Q9DCU6 
      129 AALNALQ~ state3~           26.1      26.0 state3_3   stat~         3 Q78T54 
      130 ESSSHHP~ state3~           28.3      28.2 state3_3   stat~         3 P02671 
      131 FEAHPND~ state3~           24.5      24.4 state3_3   stat~         3 G3UYD0~
      132 LPNSVLGR state3~           25.1      24.9 state3_3   stat~         3 Q8BH64 
      133 LQQDIEA~ state3~           28.9      28.7 state3_3   stat~         3 Q99M74 
      134 MELERPG~ state3~           29.1      29.0 state3_3   stat~         3 P02671 
      135 MQEVVAN~ state3~           22.2      22.1 state3_3   stat~         3 Q3U1J4 
      136 NNKDSHS~ state3~           25.3      25.2 state3_3   stat~         3 P02671 
      137 NVVIAAD~ state3~           NA        NA   state3_3   stat~         3 Q9ESL4 
      138 QGQDGLL~ state3~           25.3      25.2 state3_3   stat~         3 P57759 
      139 QHIEKAK  state3~           NA        NA   state3_3   stat~         3 A0A087~
      140 RLDETPD~ state3~           NA        NA   state3_3   stat~         3 Q61586 
      141 VEADIAG~ state3~           29.9      29.8 state3_3   stat~         3 contam~
      142 VGVTVAQ~ state3~           24.3      24.2 state3_3   stat~         3 A0A0A6~
      143 VVAGVAN~ state3~           32.5      32.3 state3_3   stat~         3 contam~
      144 YTPLYPFR state3~           23.7      23.5 state3_3   stat~         3 Q9DCU6 

# Data preparation without groups (loess normalization) for test_file_1

    Code
      prepData$SE
    Output
      # A SummarizedExperiment-tibble abstraction: 135 x 7
      # Features=15 | Samples=9 | Assays=intensity_norm, intensity
         .feature          .sample intensity_norm intensity SampleName peptide protein
         <chr>             <chr>            <dbl>     <dbl> <chr>      <chr>   <chr>  
       1 AIIEEYLHLNDMK     state1~           24.9      25.1 state1_1   AIIEEY~ A0A0J9~
       2 CLAFHDISPQAPTHFL~ state1~           28.0      27.8 state1_1   CLAFHD~ B0R1E3~
       3 EGWEYLK           state1~           24.6      24.6 state1_1   EGWEYLK Q8R404 
       4 FSLQDPPNK         state1~           NA        NA   state1_1   FSLQDP~ O70475 
       5 GPPPTDPYGRPPPYDR  state1~           NA        NA   state1_1   GPPPTD~ H3BJ30~
       6 KSQVFSTAADGQTQVE~ state1~           24.2      24.5 state1_1   KSQVFS~ P38647 
       7 LQISHEAAACITALR   state1~           25.1      25.1 state1_1   LQISHE~ A0A087~
       8 MLVDDIGDVTITNDGA~ state1~           25.3      25.2 state1_1   MLVDDI~ F2Z483~
       9 QLIVGVNK          state1~           31.3      31.3 state1_1   QLIVGV~ D3YZ68~
      10 RGEDMMHPLK        state1~           21.2      20.5 state1_1   RGEDMM~ Q9QYJ0 
      # i 125 more rows

---

    Code
      SummarizedExperiment::colData(prepData$SE)
    Output
      DataFrame with 9 rows and 1 column
                SampleName
               <character>
      state1_1    state1_1
      state1_2    state1_2
      state1_3    state1_3
      state2_1    state2_1
      state2_2    state2_2
      state2_3    state2_3
      state3_1    state3_1
      state3_2    state3_2
      state3_3    state3_3

---

    Code
      SummarizedExperiment::rowData(prepData$SE)
    Output
      DataFrame with 15 rows and 2 columns
                                        peptide                protein
                                    <character>            <character>
      AIIEEYLHLNDMK               AIIEEYLHLNDMK A0A0J9YUS5/E9PVC5/E9..
      CLAFHDISPQAPTHFLVIPK CLAFHDISPQAPTHFLVIPK          B0R1E3/P70349
      EGWEYLK                           EGWEYLK                 Q8R404
      FSLQDPPNK                       FSLQDPPNK                 O70475
      GPPPTDPYGRPPPYDR         GPPPTDPYGRPPPYDR H3BJ30/H3BJW3/H3BKW0..
      ...                                   ...                    ...
      SPLAQMEEERR                   SPLAQMEEERR   E9Q1G8/E9Q9F5/O55131
      TILTLTGVSSLEDVK           TILTLTGVSSLEDVK                 Q8CHP8
      VGEATETALTCLVEK           VGEATETALTCLVEK B1ATS4/B1ATS5/E9Q559..
      WYLTLAR                           WYLTLAR                 Q9CRA4
      YAALSDQGLDIK                 YAALSDQGLDIK                 Q9ET54

---

    Code
      SummarizedExperiment::assays(prepData$SE)$intensity
    Output
                           state1_1 state1_2 state1_3 state2_1 state2_2 state2_3
      AIIEEYLHLNDMK        25.05714 25.00897 25.00897 25.30862 25.36830 25.22092
      CLAFHDISPQAPTHFLVIPK 27.76904 27.78055 27.78055 27.76946 27.69470 27.80131
      EGWEYLK              24.58363 24.52493 24.52493 24.41524 24.66698 26.14177
      FSLQDPPNK                  NA       NA       NA       NA       NA       NA
      GPPPTDPYGRPPPYDR           NA 24.65419 24.65419       NA 24.66177       NA
      KSQVFSTAADGQTQVEIK   24.50016 24.41815 24.41815 24.41136 24.31976 24.48265
      LQISHEAAACITALR      25.10423 25.19839 25.19839 25.17251 25.26111 25.25880
      MLVDDIGDVTITNDGATILK 25.17914 24.96077 24.96077 24.99829 25.20865 24.98054
      QLIVGVNK             31.27068 31.33413 31.33413 31.22595 31.19187 31.21794
      RGEDMMHPLK           20.51265       NA       NA 20.50526 21.00883 20.50779
      SPLAQMEEERR          25.63831 25.63703 25.63703 25.49616 25.52192 25.56214
      TILTLTGVSSLEDVK            NA 22.43826 22.43826 23.23123       NA 23.03369
      VGEATETALTCLVEK      25.56135 25.57554 25.57554 25.63057 25.58600 25.49540
      WYLTLAR              23.92615       NA       NA       NA 22.70521 23.70650
      YAALSDQGLDIK         25.66946 25.84643 25.84643 25.81062 25.90850 25.85126
                           state3_1 state3_2 state3_3
      AIIEEYLHLNDMK        25.38862 25.17407 25.24943
      CLAFHDISPQAPTHFLVIPK 27.88765 28.00929 28.10054
      EGWEYLK              24.75417 24.60624 24.63681
      FSLQDPPNK                  NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA 22.93215
      KSQVFSTAADGQTQVEIK   24.38280 24.56857 24.40213
      LQISHEAAACITALR            NA 25.14672 25.11082
      MLVDDIGDVTITNDGATILK 25.11336 25.06149 25.03569
      QLIVGVNK             31.32439 31.18036 31.26893
      RGEDMMHPLK           21.51618 20.24024 21.27177
      SPLAQMEEERR          25.60741 25.59028 25.68574
      TILTLTGVSSLEDVK      23.19799       NA 22.99236
      VGEATETALTCLVEK      25.61714 25.50280 25.37474
      WYLTLAR                    NA       NA 23.64954
      YAALSDQGLDIK         25.90193 25.90039 25.89021

---

    Code
      SummarizedExperiment::assays(prepData$SE)$intensity_norm
    Output
                           state1_1 state1_2 state1_3 state2_1 state2_2 state2_3
      AIIEEYLHLNDMK        24.94314 25.13117 25.13117 25.29689 25.38346 25.06740
      CLAFHDISPQAPTHFLVIPK 28.01547 28.01547 28.01547 28.01547 28.21293 28.21293
      EGWEYLK              24.59441 24.19213 24.19213 24.59441 24.83308 26.37040
      FSLQDPPNK                  NA       NA       NA       NA       NA       NA
      GPPPTDPYGRPPPYDR           NA 24.59441 24.59441       NA 24.50869       NA
      KSQVFSTAADGQTQVEIK   24.19213 23.33253 23.33253 24.19213 24.10726 24.50869
      LQISHEAAACITALR      25.13117 25.29689 25.29689 25.13117 25.19394 25.19394
      MLVDDIGDVTITNDGATILK 25.29689 24.94314 24.94314 24.94314 25.06740 24.83308
      QLIVGVNK             31.26093 31.26093 31.26093 31.26093 31.26093 31.26093
      RGEDMMHPLK           21.15991       NA       NA 21.15991 21.15991 21.15991
      SPLAQMEEERR          25.64670 25.64670 25.64670 25.50705 25.55320 25.55320
      TILTLTGVSSLEDVK            NA 21.15991 21.15991 23.33253       NA 23.23900
      VGEATETALTCLVEK      25.50705 25.50705 25.50705 25.64670 25.70177 25.38346
      WYLTLAR              23.33253       NA       NA       NA 23.23900 24.10726
      YAALSDQGLDIK         26.18519 26.18519 26.18519 26.18519 26.37040 25.70177
                           state3_1 state3_2 state3_3
      AIIEEYLHLNDMK        25.19394 25.19394 25.45816
      CLAFHDISPQAPTHFLVIPK 27.77851 27.77851 28.38001
      EGWEYLK              24.70029 24.29398 24.99229
      FSLQDPPNK                  NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA 23.15986
      KSQVFSTAADGQTQVEIK   24.29398 23.44477 24.73094
      LQISHEAAACITALR            NA 25.02835 25.28105
      MLVDDIGDVTITNDGATILK 25.02835 24.70029 25.14082
      QLIVGVNK             31.26093 31.26093 31.26093
      RGEDMMHPLK           21.15991 21.15991 21.15991
      SPLAQMEEERR          25.42230 25.59662 25.74837
      TILTLTGVSSLEDVK      23.44477       NA 24.03545
      VGEATETALTCLVEK      25.59662 25.42230 25.58660
      WYLTLAR                    NA       NA 24.43616
      YAALSDQGLDIK         25.96294 25.96294 26.52711

---

    Code
      print(prepData$D_long, n = 200)
    Output
      # A tibble: 135 x 6
          .feature             .sample  intensity_norm intensity SampleName protein   
          <chr>                <fct>             <dbl>     <dbl> <chr>      <chr>     
        1 AIIEEYLHLNDMK        state1_1           24.9      25.1 state1_1   A0A0J9YUS~
        2 CLAFHDISPQAPTHFLVIPK state1_1           28.0      27.8 state1_1   B0R1E3/P7~
        3 EGWEYLK              state1_1           24.6      24.6 state1_1   Q8R404    
        4 FSLQDPPNK            state1_1           NA        NA   state1_1   O70475    
        5 GPPPTDPYGRPPPYDR     state1_1           NA        NA   state1_1   H3BJ30/H3~
        6 KSQVFSTAADGQTQVEIK   state1_1           24.2      24.5 state1_1   P38647    
        7 LQISHEAAACITALR      state1_1           25.1      25.1 state1_1   A0A087WPL~
        8 MLVDDIGDVTITNDGATILK state1_1           25.3      25.2 state1_1   F2Z483/P1~
        9 QLIVGVNK             state1_1           31.3      31.3 state1_1   D3YZ68/D3~
       10 RGEDMMHPLK           state1_1           21.2      20.5 state1_1   Q9QYJ0    
       11 SPLAQMEEERR          state1_1           25.6      25.6 state1_1   E9Q1G8/E9~
       12 TILTLTGVSSLEDVK      state1_1           NA        NA   state1_1   Q8CHP8    
       13 VGEATETALTCLVEK      state1_1           25.5      25.6 state1_1   B1ATS4/B1~
       14 WYLTLAR              state1_1           23.3      23.9 state1_1   Q9CRA4    
       15 YAALSDQGLDIK         state1_1           26.2      25.7 state1_1   Q9ET54    
       16 AIIEEYLHLNDMK        state1_2           25.1      25.0 state1_2   A0A0J9YUS~
       17 CLAFHDISPQAPTHFLVIPK state1_2           28.0      27.8 state1_2   B0R1E3/P7~
       18 EGWEYLK              state1_2           24.2      24.5 state1_2   Q8R404    
       19 FSLQDPPNK            state1_2           NA        NA   state1_2   O70475    
       20 GPPPTDPYGRPPPYDR     state1_2           24.6      24.7 state1_2   H3BJ30/H3~
       21 KSQVFSTAADGQTQVEIK   state1_2           23.3      24.4 state1_2   P38647    
       22 LQISHEAAACITALR      state1_2           25.3      25.2 state1_2   A0A087WPL~
       23 MLVDDIGDVTITNDGATILK state1_2           24.9      25.0 state1_2   F2Z483/P1~
       24 QLIVGVNK             state1_2           31.3      31.3 state1_2   D3YZ68/D3~
       25 RGEDMMHPLK           state1_2           NA        NA   state1_2   Q9QYJ0    
       26 SPLAQMEEERR          state1_2           25.6      25.6 state1_2   E9Q1G8/E9~
       27 TILTLTGVSSLEDVK      state1_2           21.2      22.4 state1_2   Q8CHP8    
       28 VGEATETALTCLVEK      state1_2           25.5      25.6 state1_2   B1ATS4/B1~
       29 WYLTLAR              state1_2           NA        NA   state1_2   Q9CRA4    
       30 YAALSDQGLDIK         state1_2           26.2      25.8 state1_2   Q9ET54    
       31 AIIEEYLHLNDMK        state1_3           25.1      25.0 state1_3   A0A0J9YUS~
       32 CLAFHDISPQAPTHFLVIPK state1_3           28.0      27.8 state1_3   B0R1E3/P7~
       33 EGWEYLK              state1_3           24.2      24.5 state1_3   Q8R404    
       34 FSLQDPPNK            state1_3           NA        NA   state1_3   O70475    
       35 GPPPTDPYGRPPPYDR     state1_3           24.6      24.7 state1_3   H3BJ30/H3~
       36 KSQVFSTAADGQTQVEIK   state1_3           23.3      24.4 state1_3   P38647    
       37 LQISHEAAACITALR      state1_3           25.3      25.2 state1_3   A0A087WPL~
       38 MLVDDIGDVTITNDGATILK state1_3           24.9      25.0 state1_3   F2Z483/P1~
       39 QLIVGVNK             state1_3           31.3      31.3 state1_3   D3YZ68/D3~
       40 RGEDMMHPLK           state1_3           NA        NA   state1_3   Q9QYJ0    
       41 SPLAQMEEERR          state1_3           25.6      25.6 state1_3   E9Q1G8/E9~
       42 TILTLTGVSSLEDVK      state1_3           21.2      22.4 state1_3   Q8CHP8    
       43 VGEATETALTCLVEK      state1_3           25.5      25.6 state1_3   B1ATS4/B1~
       44 WYLTLAR              state1_3           NA        NA   state1_3   Q9CRA4    
       45 YAALSDQGLDIK         state1_3           26.2      25.8 state1_3   Q9ET54    
       46 AIIEEYLHLNDMK        state2_1           25.3      25.3 state2_1   A0A0J9YUS~
       47 CLAFHDISPQAPTHFLVIPK state2_1           28.0      27.8 state2_1   B0R1E3/P7~
       48 EGWEYLK              state2_1           24.6      24.4 state2_1   Q8R404    
       49 FSLQDPPNK            state2_1           NA        NA   state2_1   O70475    
       50 GPPPTDPYGRPPPYDR     state2_1           NA        NA   state2_1   H3BJ30/H3~
       51 KSQVFSTAADGQTQVEIK   state2_1           24.2      24.4 state2_1   P38647    
       52 LQISHEAAACITALR      state2_1           25.1      25.2 state2_1   A0A087WPL~
       53 MLVDDIGDVTITNDGATILK state2_1           24.9      25.0 state2_1   F2Z483/P1~
       54 QLIVGVNK             state2_1           31.3      31.2 state2_1   D3YZ68/D3~
       55 RGEDMMHPLK           state2_1           21.2      20.5 state2_1   Q9QYJ0    
       56 SPLAQMEEERR          state2_1           25.5      25.5 state2_1   E9Q1G8/E9~
       57 TILTLTGVSSLEDVK      state2_1           23.3      23.2 state2_1   Q8CHP8    
       58 VGEATETALTCLVEK      state2_1           25.6      25.6 state2_1   B1ATS4/B1~
       59 WYLTLAR              state2_1           NA        NA   state2_1   Q9CRA4    
       60 YAALSDQGLDIK         state2_1           26.2      25.8 state2_1   Q9ET54    
       61 AIIEEYLHLNDMK        state2_2           25.4      25.4 state2_2   A0A0J9YUS~
       62 CLAFHDISPQAPTHFLVIPK state2_2           28.2      27.7 state2_2   B0R1E3/P7~
       63 EGWEYLK              state2_2           24.8      24.7 state2_2   Q8R404    
       64 FSLQDPPNK            state2_2           NA        NA   state2_2   O70475    
       65 GPPPTDPYGRPPPYDR     state2_2           24.5      24.7 state2_2   H3BJ30/H3~
       66 KSQVFSTAADGQTQVEIK   state2_2           24.1      24.3 state2_2   P38647    
       67 LQISHEAAACITALR      state2_2           25.2      25.3 state2_2   A0A087WPL~
       68 MLVDDIGDVTITNDGATILK state2_2           25.1      25.2 state2_2   F2Z483/P1~
       69 QLIVGVNK             state2_2           31.3      31.2 state2_2   D3YZ68/D3~
       70 RGEDMMHPLK           state2_2           21.2      21.0 state2_2   Q9QYJ0    
       71 SPLAQMEEERR          state2_2           25.6      25.5 state2_2   E9Q1G8/E9~
       72 TILTLTGVSSLEDVK      state2_2           NA        NA   state2_2   Q8CHP8    
       73 VGEATETALTCLVEK      state2_2           25.7      25.6 state2_2   B1ATS4/B1~
       74 WYLTLAR              state2_2           23.2      22.7 state2_2   Q9CRA4    
       75 YAALSDQGLDIK         state2_2           26.4      25.9 state2_2   Q9ET54    
       76 AIIEEYLHLNDMK        state2_3           25.1      25.2 state2_3   A0A0J9YUS~
       77 CLAFHDISPQAPTHFLVIPK state2_3           28.2      27.8 state2_3   B0R1E3/P7~
       78 EGWEYLK              state2_3           26.4      26.1 state2_3   Q8R404    
       79 FSLQDPPNK            state2_3           NA        NA   state2_3   O70475    
       80 GPPPTDPYGRPPPYDR     state2_3           NA        NA   state2_3   H3BJ30/H3~
       81 KSQVFSTAADGQTQVEIK   state2_3           24.5      24.5 state2_3   P38647    
       82 LQISHEAAACITALR      state2_3           25.2      25.3 state2_3   A0A087WPL~
       83 MLVDDIGDVTITNDGATILK state2_3           24.8      25.0 state2_3   F2Z483/P1~
       84 QLIVGVNK             state2_3           31.3      31.2 state2_3   D3YZ68/D3~
       85 RGEDMMHPLK           state2_3           21.2      20.5 state2_3   Q9QYJ0    
       86 SPLAQMEEERR          state2_3           25.6      25.6 state2_3   E9Q1G8/E9~
       87 TILTLTGVSSLEDVK      state2_3           23.2      23.0 state2_3   Q8CHP8    
       88 VGEATETALTCLVEK      state2_3           25.4      25.5 state2_3   B1ATS4/B1~
       89 WYLTLAR              state2_3           24.1      23.7 state2_3   Q9CRA4    
       90 YAALSDQGLDIK         state2_3           25.7      25.9 state2_3   Q9ET54    
       91 AIIEEYLHLNDMK        state3_1           25.2      25.4 state3_1   A0A0J9YUS~
       92 CLAFHDISPQAPTHFLVIPK state3_1           27.8      27.9 state3_1   B0R1E3/P7~
       93 EGWEYLK              state3_1           24.7      24.8 state3_1   Q8R404    
       94 FSLQDPPNK            state3_1           NA        NA   state3_1   O70475    
       95 GPPPTDPYGRPPPYDR     state3_1           NA        NA   state3_1   H3BJ30/H3~
       96 KSQVFSTAADGQTQVEIK   state3_1           24.3      24.4 state3_1   P38647    
       97 LQISHEAAACITALR      state3_1           NA        NA   state3_1   A0A087WPL~
       98 MLVDDIGDVTITNDGATILK state3_1           25.0      25.1 state3_1   F2Z483/P1~
       99 QLIVGVNK             state3_1           31.3      31.3 state3_1   D3YZ68/D3~
      100 RGEDMMHPLK           state3_1           21.2      21.5 state3_1   Q9QYJ0    
      101 SPLAQMEEERR          state3_1           25.4      25.6 state3_1   E9Q1G8/E9~
      102 TILTLTGVSSLEDVK      state3_1           23.4      23.2 state3_1   Q8CHP8    
      103 VGEATETALTCLVEK      state3_1           25.6      25.6 state3_1   B1ATS4/B1~
      104 WYLTLAR              state3_1           NA        NA   state3_1   Q9CRA4    
      105 YAALSDQGLDIK         state3_1           26.0      25.9 state3_1   Q9ET54    
      106 AIIEEYLHLNDMK        state3_2           25.2      25.2 state3_2   A0A0J9YUS~
      107 CLAFHDISPQAPTHFLVIPK state3_2           27.8      28.0 state3_2   B0R1E3/P7~
      108 EGWEYLK              state3_2           24.3      24.6 state3_2   Q8R404    
      109 FSLQDPPNK            state3_2           NA        NA   state3_2   O70475    
      110 GPPPTDPYGRPPPYDR     state3_2           NA        NA   state3_2   H3BJ30/H3~
      111 KSQVFSTAADGQTQVEIK   state3_2           23.4      24.6 state3_2   P38647    
      112 LQISHEAAACITALR      state3_2           25.0      25.1 state3_2   A0A087WPL~
      113 MLVDDIGDVTITNDGATILK state3_2           24.7      25.1 state3_2   F2Z483/P1~
      114 QLIVGVNK             state3_2           31.3      31.2 state3_2   D3YZ68/D3~
      115 RGEDMMHPLK           state3_2           21.2      20.2 state3_2   Q9QYJ0    
      116 SPLAQMEEERR          state3_2           25.6      25.6 state3_2   E9Q1G8/E9~
      117 TILTLTGVSSLEDVK      state3_2           NA        NA   state3_2   Q8CHP8    
      118 VGEATETALTCLVEK      state3_2           25.4      25.5 state3_2   B1ATS4/B1~
      119 WYLTLAR              state3_2           NA        NA   state3_2   Q9CRA4    
      120 YAALSDQGLDIK         state3_2           26.0      25.9 state3_2   Q9ET54    
      121 AIIEEYLHLNDMK        state3_3           25.5      25.2 state3_3   A0A0J9YUS~
      122 CLAFHDISPQAPTHFLVIPK state3_3           28.4      28.1 state3_3   B0R1E3/P7~
      123 EGWEYLK              state3_3           25.0      24.6 state3_3   Q8R404    
      124 FSLQDPPNK            state3_3           NA        NA   state3_3   O70475    
      125 GPPPTDPYGRPPPYDR     state3_3           23.2      22.9 state3_3   H3BJ30/H3~
      126 KSQVFSTAADGQTQVEIK   state3_3           24.7      24.4 state3_3   P38647    
      127 LQISHEAAACITALR      state3_3           25.3      25.1 state3_3   A0A087WPL~
      128 MLVDDIGDVTITNDGATILK state3_3           25.1      25.0 state3_3   F2Z483/P1~
      129 QLIVGVNK             state3_3           31.3      31.3 state3_3   D3YZ68/D3~
      130 RGEDMMHPLK           state3_3           21.2      21.3 state3_3   Q9QYJ0    
      131 SPLAQMEEERR          state3_3           25.7      25.7 state3_3   E9Q1G8/E9~
      132 TILTLTGVSSLEDVK      state3_3           24.0      23.0 state3_3   Q8CHP8    
      133 VGEATETALTCLVEK      state3_3           25.6      25.4 state3_3   B1ATS4/B1~
      134 WYLTLAR              state3_3           24.4      23.6 state3_3   Q9CRA4    
      135 YAALSDQGLDIK         state3_3           26.5      25.9 state3_3   Q9ET54    

