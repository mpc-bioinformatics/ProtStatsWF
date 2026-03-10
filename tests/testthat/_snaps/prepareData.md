# Data preparation with groups (median normalization) for test_file_1

    Code
      prepData$SE
    Output
      class: SummarizedExperiment 
      dim: 15 9 
      metadata(0):
      assays(2): intensity_norm intensity
      rownames(15): AIIEEYLHLNDMK CLAFHDISPQAPTHFLVIPK ... WYLTLAR
        YAALSDQGLDIK
      rowData names(2): peptide protein
      colnames(9): state1_1 state1_2 ... state3_2 state3_3
      colData names(3): SampleName group replicate

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
      AIIEEYLHLNDMK        25.05714       NA 25.00897 25.30862 25.36830 25.22092
      CLAFHDISPQAPTHFLVIPK 27.76904       NA 27.78055 27.76946 27.69470 27.80131
      EGWEYLK              24.58363       NA 24.52493 24.41524 24.66698 26.14177
      FSLQDPPNK                  NA       NA       NA       NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA 24.65419       NA 24.66177       NA
      KSQVFSTAADGQTQVEIK   24.50016       NA 24.41815 24.41136 24.31976 24.48265
      LQISHEAAACITALR      25.10423       NA 25.19839 25.17251 25.26111 25.25880
      MLVDDIGDVTITNDGATILK 25.17914       NA 24.96077 24.99829 25.20865 24.98054
      QLIVGVNK             31.27068       NA 31.33413 31.22595 31.19187 31.21794
      RGEDMMHPLK           20.51265       NA       NA 20.50526 21.00883 20.50779
      SPLAQMEEERR          25.63831       NA 25.63703 25.49616 25.52192 25.56214
      TILTLTGVSSLEDVK            NA       NA 22.43826 23.23123       NA 23.03369
      VGEATETALTCLVEK      25.56135       NA 25.57554 25.63057 25.58600 25.49540
      WYLTLAR              23.92615       NA       NA       NA 22.70521 23.70650
      YAALSDQGLDIK         25.66946       NA 25.84643 25.81062 25.90850 25.85126
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
      AIIEEYLHLNDMK              NA       NA       NA       NA       NA       NA
      CLAFHDISPQAPTHFLVIPK       NA       NA       NA       NA       NA       NA
      EGWEYLK                    NA       NA       NA       NA       NA       NA
      FSLQDPPNK                  NA       NA       NA       NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA       NA       NA       NA       NA
      KSQVFSTAADGQTQVEIK         NA       NA       NA       NA       NA       NA
      LQISHEAAACITALR            NA       NA       NA       NA       NA       NA
      MLVDDIGDVTITNDGATILK       NA       NA       NA       NA       NA       NA
      QLIVGVNK                   NA       NA       NA       NA       NA       NA
      RGEDMMHPLK                 NA       NA       NA       NA       NA       NA
      SPLAQMEEERR                NA       NA       NA       NA       NA       NA
      TILTLTGVSSLEDVK            NA       NA       NA       NA       NA       NA
      VGEATETALTCLVEK            NA       NA       NA       NA       NA       NA
      WYLTLAR                    NA       NA       NA       NA       NA       NA
      YAALSDQGLDIK               NA       NA       NA       NA       NA       NA
                           state3_1 state3_2 state3_3
      AIIEEYLHLNDMK              NA       NA       NA
      CLAFHDISPQAPTHFLVIPK       NA       NA       NA
      EGWEYLK                    NA       NA       NA
      FSLQDPPNK                  NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA       NA
      KSQVFSTAADGQTQVEIK         NA       NA       NA
      LQISHEAAACITALR            NA       NA       NA
      MLVDDIGDVTITNDGATILK       NA       NA       NA
      QLIVGVNK                   NA       NA       NA
      RGEDMMHPLK                 NA       NA       NA
      SPLAQMEEERR                NA       NA       NA
      TILTLTGVSSLEDVK            NA       NA       NA
      VGEATETALTCLVEK            NA       NA       NA
      WYLTLAR                    NA       NA       NA
      YAALSDQGLDIK               NA       NA       NA

---

    Code
      print(prepData$D_long, n = 200)
    Output
      # A tibble: 135 x 8
          .feature .sample intensity_norm intensity SampleName group replicate protein
          <chr>    <fct>            <dbl>     <dbl> <chr>      <chr>     <dbl> <chr>  
        1 AIIEEYL~ state1~             NA      25.1 state1_1   stat~         1 A0A0J9~
        2 CLAFHDI~ state1~             NA      27.8 state1_1   stat~         1 B0R1E3~
        3 EGWEYLK  state1~             NA      24.6 state1_1   stat~         1 Q8R404 
        4 FSLQDPP~ state1~             NA      NA   state1_1   stat~         1 O70475 
        5 GPPPTDP~ state1~             NA      NA   state1_1   stat~         1 H3BJ30~
        6 KSQVFST~ state1~             NA      24.5 state1_1   stat~         1 P38647 
        7 LQISHEA~ state1~             NA      25.1 state1_1   stat~         1 A0A087~
        8 MLVDDIG~ state1~             NA      25.2 state1_1   stat~         1 F2Z483~
        9 QLIVGVNK state1~             NA      31.3 state1_1   stat~         1 D3YZ68~
       10 RGEDMMH~ state1~             NA      20.5 state1_1   stat~         1 Q9QYJ0 
       11 SPLAQME~ state1~             NA      25.6 state1_1   stat~         1 E9Q1G8~
       12 TILTLTG~ state1~             NA      NA   state1_1   stat~         1 Q8CHP8 
       13 VGEATET~ state1~             NA      25.6 state1_1   stat~         1 B1ATS4~
       14 WYLTLAR  state1~             NA      23.9 state1_1   stat~         1 Q9CRA4 
       15 YAALSDQ~ state1~             NA      25.7 state1_1   stat~         1 Q9ET54 
       16 AIIEEYL~ state1~             NA      NA   state1_2   stat~         2 A0A0J9~
       17 CLAFHDI~ state1~             NA      NA   state1_2   stat~         2 B0R1E3~
       18 EGWEYLK  state1~             NA      NA   state1_2   stat~         2 Q8R404 
       19 FSLQDPP~ state1~             NA      NA   state1_2   stat~         2 O70475 
       20 GPPPTDP~ state1~             NA      NA   state1_2   stat~         2 H3BJ30~
       21 KSQVFST~ state1~             NA      NA   state1_2   stat~         2 P38647 
       22 LQISHEA~ state1~             NA      NA   state1_2   stat~         2 A0A087~
       23 MLVDDIG~ state1~             NA      NA   state1_2   stat~         2 F2Z483~
       24 QLIVGVNK state1~             NA      NA   state1_2   stat~         2 D3YZ68~
       25 RGEDMMH~ state1~             NA      NA   state1_2   stat~         2 Q9QYJ0 
       26 SPLAQME~ state1~             NA      NA   state1_2   stat~         2 E9Q1G8~
       27 TILTLTG~ state1~             NA      NA   state1_2   stat~         2 Q8CHP8 
       28 VGEATET~ state1~             NA      NA   state1_2   stat~         2 B1ATS4~
       29 WYLTLAR  state1~             NA      NA   state1_2   stat~         2 Q9CRA4 
       30 YAALSDQ~ state1~             NA      NA   state1_2   stat~         2 Q9ET54 
       31 AIIEEYL~ state1~             NA      25.0 state1_3   stat~         3 A0A0J9~
       32 CLAFHDI~ state1~             NA      27.8 state1_3   stat~         3 B0R1E3~
       33 EGWEYLK  state1~             NA      24.5 state1_3   stat~         3 Q8R404 
       34 FSLQDPP~ state1~             NA      NA   state1_3   stat~         3 O70475 
       35 GPPPTDP~ state1~             NA      24.7 state1_3   stat~         3 H3BJ30~
       36 KSQVFST~ state1~             NA      24.4 state1_3   stat~         3 P38647 
       37 LQISHEA~ state1~             NA      25.2 state1_3   stat~         3 A0A087~
       38 MLVDDIG~ state1~             NA      25.0 state1_3   stat~         3 F2Z483~
       39 QLIVGVNK state1~             NA      31.3 state1_3   stat~         3 D3YZ68~
       40 RGEDMMH~ state1~             NA      NA   state1_3   stat~         3 Q9QYJ0 
       41 SPLAQME~ state1~             NA      25.6 state1_3   stat~         3 E9Q1G8~
       42 TILTLTG~ state1~             NA      22.4 state1_3   stat~         3 Q8CHP8 
       43 VGEATET~ state1~             NA      25.6 state1_3   stat~         3 B1ATS4~
       44 WYLTLAR  state1~             NA      NA   state1_3   stat~         3 Q9CRA4 
       45 YAALSDQ~ state1~             NA      25.8 state1_3   stat~         3 Q9ET54 
       46 AIIEEYL~ state2~             NA      25.3 state2_1   stat~         1 A0A0J9~
       47 CLAFHDI~ state2~             NA      27.8 state2_1   stat~         1 B0R1E3~
       48 EGWEYLK  state2~             NA      24.4 state2_1   stat~         1 Q8R404 
       49 FSLQDPP~ state2~             NA      NA   state2_1   stat~         1 O70475 
       50 GPPPTDP~ state2~             NA      NA   state2_1   stat~         1 H3BJ30~
       51 KSQVFST~ state2~             NA      24.4 state2_1   stat~         1 P38647 
       52 LQISHEA~ state2~             NA      25.2 state2_1   stat~         1 A0A087~
       53 MLVDDIG~ state2~             NA      25.0 state2_1   stat~         1 F2Z483~
       54 QLIVGVNK state2~             NA      31.2 state2_1   stat~         1 D3YZ68~
       55 RGEDMMH~ state2~             NA      20.5 state2_1   stat~         1 Q9QYJ0 
       56 SPLAQME~ state2~             NA      25.5 state2_1   stat~         1 E9Q1G8~
       57 TILTLTG~ state2~             NA      23.2 state2_1   stat~         1 Q8CHP8 
       58 VGEATET~ state2~             NA      25.6 state2_1   stat~         1 B1ATS4~
       59 WYLTLAR  state2~             NA      NA   state2_1   stat~         1 Q9CRA4 
       60 YAALSDQ~ state2~             NA      25.8 state2_1   stat~         1 Q9ET54 
       61 AIIEEYL~ state2~             NA      25.4 state2_2   stat~         2 A0A0J9~
       62 CLAFHDI~ state2~             NA      27.7 state2_2   stat~         2 B0R1E3~
       63 EGWEYLK  state2~             NA      24.7 state2_2   stat~         2 Q8R404 
       64 FSLQDPP~ state2~             NA      NA   state2_2   stat~         2 O70475 
       65 GPPPTDP~ state2~             NA      24.7 state2_2   stat~         2 H3BJ30~
       66 KSQVFST~ state2~             NA      24.3 state2_2   stat~         2 P38647 
       67 LQISHEA~ state2~             NA      25.3 state2_2   stat~         2 A0A087~
       68 MLVDDIG~ state2~             NA      25.2 state2_2   stat~         2 F2Z483~
       69 QLIVGVNK state2~             NA      31.2 state2_2   stat~         2 D3YZ68~
       70 RGEDMMH~ state2~             NA      21.0 state2_2   stat~         2 Q9QYJ0 
       71 SPLAQME~ state2~             NA      25.5 state2_2   stat~         2 E9Q1G8~
       72 TILTLTG~ state2~             NA      NA   state2_2   stat~         2 Q8CHP8 
       73 VGEATET~ state2~             NA      25.6 state2_2   stat~         2 B1ATS4~
       74 WYLTLAR  state2~             NA      22.7 state2_2   stat~         2 Q9CRA4 
       75 YAALSDQ~ state2~             NA      25.9 state2_2   stat~         2 Q9ET54 
       76 AIIEEYL~ state2~             NA      25.2 state2_3   stat~         3 A0A0J9~
       77 CLAFHDI~ state2~             NA      27.8 state2_3   stat~         3 B0R1E3~
       78 EGWEYLK  state2~             NA      26.1 state2_3   stat~         3 Q8R404 
       79 FSLQDPP~ state2~             NA      NA   state2_3   stat~         3 O70475 
       80 GPPPTDP~ state2~             NA      NA   state2_3   stat~         3 H3BJ30~
       81 KSQVFST~ state2~             NA      24.5 state2_3   stat~         3 P38647 
       82 LQISHEA~ state2~             NA      25.3 state2_3   stat~         3 A0A087~
       83 MLVDDIG~ state2~             NA      25.0 state2_3   stat~         3 F2Z483~
       84 QLIVGVNK state2~             NA      31.2 state2_3   stat~         3 D3YZ68~
       85 RGEDMMH~ state2~             NA      20.5 state2_3   stat~         3 Q9QYJ0 
       86 SPLAQME~ state2~             NA      25.6 state2_3   stat~         3 E9Q1G8~
       87 TILTLTG~ state2~             NA      23.0 state2_3   stat~         3 Q8CHP8 
       88 VGEATET~ state2~             NA      25.5 state2_3   stat~         3 B1ATS4~
       89 WYLTLAR  state2~             NA      23.7 state2_3   stat~         3 Q9CRA4 
       90 YAALSDQ~ state2~             NA      25.9 state2_3   stat~         3 Q9ET54 
       91 AIIEEYL~ state3~             NA      25.4 state3_1   stat~         1 A0A0J9~
       92 CLAFHDI~ state3~             NA      27.9 state3_1   stat~         1 B0R1E3~
       93 EGWEYLK  state3~             NA      24.8 state3_1   stat~         1 Q8R404 
       94 FSLQDPP~ state3~             NA      NA   state3_1   stat~         1 O70475 
       95 GPPPTDP~ state3~             NA      NA   state3_1   stat~         1 H3BJ30~
       96 KSQVFST~ state3~             NA      24.4 state3_1   stat~         1 P38647 
       97 LQISHEA~ state3~             NA      NA   state3_1   stat~         1 A0A087~
       98 MLVDDIG~ state3~             NA      25.1 state3_1   stat~         1 F2Z483~
       99 QLIVGVNK state3~             NA      31.3 state3_1   stat~         1 D3YZ68~
      100 RGEDMMH~ state3~             NA      21.5 state3_1   stat~         1 Q9QYJ0 
      101 SPLAQME~ state3~             NA      25.6 state3_1   stat~         1 E9Q1G8~
      102 TILTLTG~ state3~             NA      23.2 state3_1   stat~         1 Q8CHP8 
      103 VGEATET~ state3~             NA      25.6 state3_1   stat~         1 B1ATS4~
      104 WYLTLAR  state3~             NA      NA   state3_1   stat~         1 Q9CRA4 
      105 YAALSDQ~ state3~             NA      25.9 state3_1   stat~         1 Q9ET54 
      106 AIIEEYL~ state3~             NA      25.2 state3_2   stat~         2 A0A0J9~
      107 CLAFHDI~ state3~             NA      28.0 state3_2   stat~         2 B0R1E3~
      108 EGWEYLK  state3~             NA      24.6 state3_2   stat~         2 Q8R404 
      109 FSLQDPP~ state3~             NA      NA   state3_2   stat~         2 O70475 
      110 GPPPTDP~ state3~             NA      NA   state3_2   stat~         2 H3BJ30~
      111 KSQVFST~ state3~             NA      24.6 state3_2   stat~         2 P38647 
      112 LQISHEA~ state3~             NA      25.1 state3_2   stat~         2 A0A087~
      113 MLVDDIG~ state3~             NA      25.1 state3_2   stat~         2 F2Z483~
      114 QLIVGVNK state3~             NA      31.2 state3_2   stat~         2 D3YZ68~
      115 RGEDMMH~ state3~             NA      20.2 state3_2   stat~         2 Q9QYJ0 
      116 SPLAQME~ state3~             NA      25.6 state3_2   stat~         2 E9Q1G8~
      117 TILTLTG~ state3~             NA      NA   state3_2   stat~         2 Q8CHP8 
      118 VGEATET~ state3~             NA      25.5 state3_2   stat~         2 B1ATS4~
      119 WYLTLAR  state3~             NA      NA   state3_2   stat~         2 Q9CRA4 
      120 YAALSDQ~ state3~             NA      25.9 state3_2   stat~         2 Q9ET54 
      121 AIIEEYL~ state3~             NA      25.2 state3_3   stat~         3 A0A0J9~
      122 CLAFHDI~ state3~             NA      28.1 state3_3   stat~         3 B0R1E3~
      123 EGWEYLK  state3~             NA      24.6 state3_3   stat~         3 Q8R404 
      124 FSLQDPP~ state3~             NA      NA   state3_3   stat~         3 O70475 
      125 GPPPTDP~ state3~             NA      22.9 state3_3   stat~         3 H3BJ30~
      126 KSQVFST~ state3~             NA      24.4 state3_3   stat~         3 P38647 
      127 LQISHEA~ state3~             NA      25.1 state3_3   stat~         3 A0A087~
      128 MLVDDIG~ state3~             NA      25.0 state3_3   stat~         3 F2Z483~
      129 QLIVGVNK state3~             NA      31.3 state3_3   stat~         3 D3YZ68~
      130 RGEDMMH~ state3~             NA      21.3 state3_3   stat~         3 Q9QYJ0 
      131 SPLAQME~ state3~             NA      25.7 state3_3   stat~         3 E9Q1G8~
      132 TILTLTG~ state3~             NA      23.0 state3_3   stat~         3 Q8CHP8 
      133 VGEATET~ state3~             NA      25.4 state3_3   stat~         3 B1ATS4~
      134 WYLTLAR  state3~             NA      23.6 state3_3   stat~         3 Q9CRA4 
      135 YAALSDQ~ state3~             NA      25.9 state3_3   stat~         3 Q9ET54 

# Data preparation without groups (loess normalization) for test_file_1

    Code
      prepData$SE
    Output
      class: SummarizedExperiment 
      dim: 15 9 
      metadata(0):
      assays(2): intensity_norm intensity
      rownames(15): AIIEEYLHLNDMK CLAFHDISPQAPTHFLVIPK ... WYLTLAR
        YAALSDQGLDIK
      rowData names(2): peptide protein
      colnames(9): state1_1 state1_2 ... state3_2 state3_3
      colData names(1): SampleName

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
      AIIEEYLHLNDMK        25.05714       NA 25.00897 25.30862 25.36830 25.22092
      CLAFHDISPQAPTHFLVIPK 27.76904       NA 27.78055 27.76946 27.69470 27.80131
      EGWEYLK              24.58363       NA 24.52493 24.41524 24.66698 26.14177
      FSLQDPPNK                  NA       NA       NA       NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA 24.65419       NA 24.66177       NA
      KSQVFSTAADGQTQVEIK   24.50016       NA 24.41815 24.41136 24.31976 24.48265
      LQISHEAAACITALR      25.10423       NA 25.19839 25.17251 25.26111 25.25880
      MLVDDIGDVTITNDGATILK 25.17914       NA 24.96077 24.99829 25.20865 24.98054
      QLIVGVNK             31.27068       NA 31.33413 31.22595 31.19187 31.21794
      RGEDMMHPLK           20.51265       NA       NA 20.50526 21.00883 20.50779
      SPLAQMEEERR          25.63831       NA 25.63703 25.49616 25.52192 25.56214
      TILTLTGVSSLEDVK            NA       NA 22.43826 23.23123       NA 23.03369
      VGEATETALTCLVEK      25.56135       NA 25.57554 25.63057 25.58600 25.49540
      WYLTLAR              23.92615       NA       NA       NA 22.70521 23.70650
      YAALSDQGLDIK         25.66946       NA 25.84643 25.81062 25.90850 25.85126
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
      AIIEEYLHLNDMK        25.17401       NA 25.13042 25.39210 25.30028 25.24947
      CLAFHDISPQAPTHFLVIPK 27.85402       NA 27.85277 27.85074 27.84997 27.84925
      EGWEYLK              24.71794       NA 24.69667 24.49887 24.64731 26.18813
      FSLQDPPNK                  NA       NA       NA       NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA 24.88293       NA 25.10851       NA
      KSQVFSTAADGQTQVEIK   24.53713       NA 24.52131 24.49102 24.49592 24.45302
      LQISHEAAACITALR      25.24380       NA 25.34360 25.26213 25.19427 25.28867
      MLVDDIGDVTITNDGATILK 25.38873       NA 25.16839 25.09942 25.14052 25.01702
      QLIVGVNK             31.25218       NA 31.25251 31.25293 31.25321 31.25333
      RGEDMMHPLK           20.79093       NA       NA 20.78549 20.79299 20.79050
      SPLAQMEEERR          25.68150       NA 25.61966 25.53219 25.44982 25.56617
      TILTLTGVSSLEDVK            NA       NA 22.98156 22.99282       NA 22.96825
      VGEATETALTCLVEK      25.60970       NA 25.57355 25.67214 25.51303 25.50290
      WYLTLAR              23.58371       NA       NA       NA 23.58957 23.63101
      YAALSDQGLDIK         25.69314       NA 25.75582 25.82146 25.84741 25.83440
                           state3_1 state3_2 state3_3
      AIIEEYLHLNDMK        25.37418 25.25304 25.30463
      CLAFHDISPQAPTHFLVIPK 27.85137 27.85151 27.85021
      EGWEYLK              24.80531 24.60544 24.77626
      FSLQDPPNK                  NA       NA       NA
      GPPPTDPYGRPPPYDR           NA       NA 22.90047
      KSQVFSTAADGQTQVEIK   24.48873 24.48872 24.53020
      LQISHEAAACITALR            NA 25.23096 25.17657
      MLVDDIGDVTITNDGATILK 25.12461 25.15537 25.13651
      QLIVGVNK             31.25287 31.25272 31.25243
      RGEDMMHPLK           20.78885 20.79058 20.79478
      SPLAQMEEERR          25.54611 25.57782 25.64942
      TILTLTGVSSLEDVK      22.98763       NA 22.96895
      VGEATETALTCLVEK      25.56196 25.50430 25.35058
      WYLTLAR                    NA       NA 23.62305
      YAALSDQGLDIK         25.81049 25.81518 25.79201

---

    Code
      print(prepData$D_long, n = 200)
    Output
      # A tibble: 135 x 6
          .feature             .sample  intensity_norm intensity SampleName protein   
          <chr>                <fct>             <dbl>     <dbl> <chr>      <chr>     
        1 AIIEEYLHLNDMK        state1_1           25.2      25.1 state1_1   A0A0J9YUS~
        2 CLAFHDISPQAPTHFLVIPK state1_1           27.9      27.8 state1_1   B0R1E3/P7~
        3 EGWEYLK              state1_1           24.7      24.6 state1_1   Q8R404    
        4 FSLQDPPNK            state1_1           NA        NA   state1_1   O70475    
        5 GPPPTDPYGRPPPYDR     state1_1           NA        NA   state1_1   H3BJ30/H3~
        6 KSQVFSTAADGQTQVEIK   state1_1           24.5      24.5 state1_1   P38647    
        7 LQISHEAAACITALR      state1_1           25.2      25.1 state1_1   A0A087WPL~
        8 MLVDDIGDVTITNDGATILK state1_1           25.4      25.2 state1_1   F2Z483/P1~
        9 QLIVGVNK             state1_1           31.3      31.3 state1_1   D3YZ68/D3~
       10 RGEDMMHPLK           state1_1           20.8      20.5 state1_1   Q9QYJ0    
       11 SPLAQMEEERR          state1_1           25.7      25.6 state1_1   E9Q1G8/E9~
       12 TILTLTGVSSLEDVK      state1_1           NA        NA   state1_1   Q8CHP8    
       13 VGEATETALTCLVEK      state1_1           25.6      25.6 state1_1   B1ATS4/B1~
       14 WYLTLAR              state1_1           23.6      23.9 state1_1   Q9CRA4    
       15 YAALSDQGLDIK         state1_1           25.7      25.7 state1_1   Q9ET54    
       16 AIIEEYLHLNDMK        state1_2           NA        NA   state1_2   A0A0J9YUS~
       17 CLAFHDISPQAPTHFLVIPK state1_2           NA        NA   state1_2   B0R1E3/P7~
       18 EGWEYLK              state1_2           NA        NA   state1_2   Q8R404    
       19 FSLQDPPNK            state1_2           NA        NA   state1_2   O70475    
       20 GPPPTDPYGRPPPYDR     state1_2           NA        NA   state1_2   H3BJ30/H3~
       21 KSQVFSTAADGQTQVEIK   state1_2           NA        NA   state1_2   P38647    
       22 LQISHEAAACITALR      state1_2           NA        NA   state1_2   A0A087WPL~
       23 MLVDDIGDVTITNDGATILK state1_2           NA        NA   state1_2   F2Z483/P1~
       24 QLIVGVNK             state1_2           NA        NA   state1_2   D3YZ68/D3~
       25 RGEDMMHPLK           state1_2           NA        NA   state1_2   Q9QYJ0    
       26 SPLAQMEEERR          state1_2           NA        NA   state1_2   E9Q1G8/E9~
       27 TILTLTGVSSLEDVK      state1_2           NA        NA   state1_2   Q8CHP8    
       28 VGEATETALTCLVEK      state1_2           NA        NA   state1_2   B1ATS4/B1~
       29 WYLTLAR              state1_2           NA        NA   state1_2   Q9CRA4    
       30 YAALSDQGLDIK         state1_2           NA        NA   state1_2   Q9ET54    
       31 AIIEEYLHLNDMK        state1_3           25.1      25.0 state1_3   A0A0J9YUS~
       32 CLAFHDISPQAPTHFLVIPK state1_3           27.9      27.8 state1_3   B0R1E3/P7~
       33 EGWEYLK              state1_3           24.7      24.5 state1_3   Q8R404    
       34 FSLQDPPNK            state1_3           NA        NA   state1_3   O70475    
       35 GPPPTDPYGRPPPYDR     state1_3           24.9      24.7 state1_3   H3BJ30/H3~
       36 KSQVFSTAADGQTQVEIK   state1_3           24.5      24.4 state1_3   P38647    
       37 LQISHEAAACITALR      state1_3           25.3      25.2 state1_3   A0A087WPL~
       38 MLVDDIGDVTITNDGATILK state1_3           25.2      25.0 state1_3   F2Z483/P1~
       39 QLIVGVNK             state1_3           31.3      31.3 state1_3   D3YZ68/D3~
       40 RGEDMMHPLK           state1_3           NA        NA   state1_3   Q9QYJ0    
       41 SPLAQMEEERR          state1_3           25.6      25.6 state1_3   E9Q1G8/E9~
       42 TILTLTGVSSLEDVK      state1_3           23.0      22.4 state1_3   Q8CHP8    
       43 VGEATETALTCLVEK      state1_3           25.6      25.6 state1_3   B1ATS4/B1~
       44 WYLTLAR              state1_3           NA        NA   state1_3   Q9CRA4    
       45 YAALSDQGLDIK         state1_3           25.8      25.8 state1_3   Q9ET54    
       46 AIIEEYLHLNDMK        state2_1           25.4      25.3 state2_1   A0A0J9YUS~
       47 CLAFHDISPQAPTHFLVIPK state2_1           27.9      27.8 state2_1   B0R1E3/P7~
       48 EGWEYLK              state2_1           24.5      24.4 state2_1   Q8R404    
       49 FSLQDPPNK            state2_1           NA        NA   state2_1   O70475    
       50 GPPPTDPYGRPPPYDR     state2_1           NA        NA   state2_1   H3BJ30/H3~
       51 KSQVFSTAADGQTQVEIK   state2_1           24.5      24.4 state2_1   P38647    
       52 LQISHEAAACITALR      state2_1           25.3      25.2 state2_1   A0A087WPL~
       53 MLVDDIGDVTITNDGATILK state2_1           25.1      25.0 state2_1   F2Z483/P1~
       54 QLIVGVNK             state2_1           31.3      31.2 state2_1   D3YZ68/D3~
       55 RGEDMMHPLK           state2_1           20.8      20.5 state2_1   Q9QYJ0    
       56 SPLAQMEEERR          state2_1           25.5      25.5 state2_1   E9Q1G8/E9~
       57 TILTLTGVSSLEDVK      state2_1           23.0      23.2 state2_1   Q8CHP8    
       58 VGEATETALTCLVEK      state2_1           25.7      25.6 state2_1   B1ATS4/B1~
       59 WYLTLAR              state2_1           NA        NA   state2_1   Q9CRA4    
       60 YAALSDQGLDIK         state2_1           25.8      25.8 state2_1   Q9ET54    
       61 AIIEEYLHLNDMK        state2_2           25.3      25.4 state2_2   A0A0J9YUS~
       62 CLAFHDISPQAPTHFLVIPK state2_2           27.8      27.7 state2_2   B0R1E3/P7~
       63 EGWEYLK              state2_2           24.6      24.7 state2_2   Q8R404    
       64 FSLQDPPNK            state2_2           NA        NA   state2_2   O70475    
       65 GPPPTDPYGRPPPYDR     state2_2           25.1      24.7 state2_2   H3BJ30/H3~
       66 KSQVFSTAADGQTQVEIK   state2_2           24.5      24.3 state2_2   P38647    
       67 LQISHEAAACITALR      state2_2           25.2      25.3 state2_2   A0A087WPL~
       68 MLVDDIGDVTITNDGATILK state2_2           25.1      25.2 state2_2   F2Z483/P1~
       69 QLIVGVNK             state2_2           31.3      31.2 state2_2   D3YZ68/D3~
       70 RGEDMMHPLK           state2_2           20.8      21.0 state2_2   Q9QYJ0    
       71 SPLAQMEEERR          state2_2           25.4      25.5 state2_2   E9Q1G8/E9~
       72 TILTLTGVSSLEDVK      state2_2           NA        NA   state2_2   Q8CHP8    
       73 VGEATETALTCLVEK      state2_2           25.5      25.6 state2_2   B1ATS4/B1~
       74 WYLTLAR              state2_2           23.6      22.7 state2_2   Q9CRA4    
       75 YAALSDQGLDIK         state2_2           25.8      25.9 state2_2   Q9ET54    
       76 AIIEEYLHLNDMK        state2_3           25.2      25.2 state2_3   A0A0J9YUS~
       77 CLAFHDISPQAPTHFLVIPK state2_3           27.8      27.8 state2_3   B0R1E3/P7~
       78 EGWEYLK              state2_3           26.2      26.1 state2_3   Q8R404    
       79 FSLQDPPNK            state2_3           NA        NA   state2_3   O70475    
       80 GPPPTDPYGRPPPYDR     state2_3           NA        NA   state2_3   H3BJ30/H3~
       81 KSQVFSTAADGQTQVEIK   state2_3           24.5      24.5 state2_3   P38647    
       82 LQISHEAAACITALR      state2_3           25.3      25.3 state2_3   A0A087WPL~
       83 MLVDDIGDVTITNDGATILK state2_3           25.0      25.0 state2_3   F2Z483/P1~
       84 QLIVGVNK             state2_3           31.3      31.2 state2_3   D3YZ68/D3~
       85 RGEDMMHPLK           state2_3           20.8      20.5 state2_3   Q9QYJ0    
       86 SPLAQMEEERR          state2_3           25.6      25.6 state2_3   E9Q1G8/E9~
       87 TILTLTGVSSLEDVK      state2_3           23.0      23.0 state2_3   Q8CHP8    
       88 VGEATETALTCLVEK      state2_3           25.5      25.5 state2_3   B1ATS4/B1~
       89 WYLTLAR              state2_3           23.6      23.7 state2_3   Q9CRA4    
       90 YAALSDQGLDIK         state2_3           25.8      25.9 state2_3   Q9ET54    
       91 AIIEEYLHLNDMK        state3_1           25.4      25.4 state3_1   A0A0J9YUS~
       92 CLAFHDISPQAPTHFLVIPK state3_1           27.9      27.9 state3_1   B0R1E3/P7~
       93 EGWEYLK              state3_1           24.8      24.8 state3_1   Q8R404    
       94 FSLQDPPNK            state3_1           NA        NA   state3_1   O70475    
       95 GPPPTDPYGRPPPYDR     state3_1           NA        NA   state3_1   H3BJ30/H3~
       96 KSQVFSTAADGQTQVEIK   state3_1           24.5      24.4 state3_1   P38647    
       97 LQISHEAAACITALR      state3_1           NA        NA   state3_1   A0A087WPL~
       98 MLVDDIGDVTITNDGATILK state3_1           25.1      25.1 state3_1   F2Z483/P1~
       99 QLIVGVNK             state3_1           31.3      31.3 state3_1   D3YZ68/D3~
      100 RGEDMMHPLK           state3_1           20.8      21.5 state3_1   Q9QYJ0    
      101 SPLAQMEEERR          state3_1           25.5      25.6 state3_1   E9Q1G8/E9~
      102 TILTLTGVSSLEDVK      state3_1           23.0      23.2 state3_1   Q8CHP8    
      103 VGEATETALTCLVEK      state3_1           25.6      25.6 state3_1   B1ATS4/B1~
      104 WYLTLAR              state3_1           NA        NA   state3_1   Q9CRA4    
      105 YAALSDQGLDIK         state3_1           25.8      25.9 state3_1   Q9ET54    
      106 AIIEEYLHLNDMK        state3_2           25.3      25.2 state3_2   A0A0J9YUS~
      107 CLAFHDISPQAPTHFLVIPK state3_2           27.9      28.0 state3_2   B0R1E3/P7~
      108 EGWEYLK              state3_2           24.6      24.6 state3_2   Q8R404    
      109 FSLQDPPNK            state3_2           NA        NA   state3_2   O70475    
      110 GPPPTDPYGRPPPYDR     state3_2           NA        NA   state3_2   H3BJ30/H3~
      111 KSQVFSTAADGQTQVEIK   state3_2           24.5      24.6 state3_2   P38647    
      112 LQISHEAAACITALR      state3_2           25.2      25.1 state3_2   A0A087WPL~
      113 MLVDDIGDVTITNDGATILK state3_2           25.2      25.1 state3_2   F2Z483/P1~
      114 QLIVGVNK             state3_2           31.3      31.2 state3_2   D3YZ68/D3~
      115 RGEDMMHPLK           state3_2           20.8      20.2 state3_2   Q9QYJ0    
      116 SPLAQMEEERR          state3_2           25.6      25.6 state3_2   E9Q1G8/E9~
      117 TILTLTGVSSLEDVK      state3_2           NA        NA   state3_2   Q8CHP8    
      118 VGEATETALTCLVEK      state3_2           25.5      25.5 state3_2   B1ATS4/B1~
      119 WYLTLAR              state3_2           NA        NA   state3_2   Q9CRA4    
      120 YAALSDQGLDIK         state3_2           25.8      25.9 state3_2   Q9ET54    
      121 AIIEEYLHLNDMK        state3_3           25.3      25.2 state3_3   A0A0J9YUS~
      122 CLAFHDISPQAPTHFLVIPK state3_3           27.9      28.1 state3_3   B0R1E3/P7~
      123 EGWEYLK              state3_3           24.8      24.6 state3_3   Q8R404    
      124 FSLQDPPNK            state3_3           NA        NA   state3_3   O70475    
      125 GPPPTDPYGRPPPYDR     state3_3           22.9      22.9 state3_3   H3BJ30/H3~
      126 KSQVFSTAADGQTQVEIK   state3_3           24.5      24.4 state3_3   P38647    
      127 LQISHEAAACITALR      state3_3           25.2      25.1 state3_3   A0A087WPL~
      128 MLVDDIGDVTITNDGATILK state3_3           25.1      25.0 state3_3   F2Z483/P1~
      129 QLIVGVNK             state3_3           31.3      31.3 state3_3   D3YZ68/D3~
      130 RGEDMMHPLK           state3_3           20.8      21.3 state3_3   Q9QYJ0    
      131 SPLAQMEEERR          state3_3           25.6      25.7 state3_3   E9Q1G8/E9~
      132 TILTLTGVSSLEDVK      state3_3           23.0      23.0 state3_3   Q8CHP8    
      133 VGEATETALTCLVEK      state3_3           25.4      25.4 state3_3   B1ATS4/B1~
      134 WYLTLAR              state3_3           23.6      23.6 state3_3   Q9CRA4    
      135 YAALSDQGLDIK         state3_3           25.8      25.9 state3_3   Q9ET54    

