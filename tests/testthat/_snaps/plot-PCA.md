# Test standard PCA plot

    Code
      PCA$D_PCA_plot
    Output
                     PCx         PCy  label Sample PatientID Group Gender Age
      HCC_1  -24.2260058  -5.8796975  HCC_1  HCC_1         1   HCC Female  57
      HCC_2   -0.2701998  11.7565454  HCC_2  HCC_2         2   HCC   Male  71
      HCC_3  -32.8285745   2.9649821  HCC_3  HCC_3         3   HCC Female  42
      HCC_4   -8.0581633   3.2113162  HCC_4  HCC_4         4   HCC   Male  73
      HCC_5    5.4800309   4.1258811  HCC_5  HCC_5         5   HCC   Male  51
      HCC_6   -9.0077163  -0.3124016  HCC_6  HCC_6         6   HCC   Male  21
      HCC_7  -16.6865955  -0.9581829  HCC_7  HCC_7         7   HCC   Male  71
      HCC_8   -0.7120405   7.0650912  HCC_8  HCC_8         8   HCC   Male  74
      HCC_9   -8.8132365   4.0503016  HCC_9  HCC_9         9   HCC   Male  65
      HCC_10 -46.4426474 -16.5206031 HCC_10 HCC_10        10   HCC Female  59
      HCC_11 -13.7105749   2.8033642 HCC_11 HCC_11        11   HCC Female  76
      HCC_12 -14.0413060   0.3658427 HCC_12 HCC_12        12   HCC   Male  53
      HCC_13 -14.5881346   3.4584455 HCC_13 HCC_13        13   HCC   Male  62
      HCC_14   3.3520552   6.3015778 HCC_14 HCC_14        14   HCC   Male  72
      HCC_15 -19.2159911   2.4651543 HCC_15 HCC_15        15   HCC Female  67
      HCC_16 -10.6083777   2.7896151 HCC_16 HCC_16        16   HCC   Male  75
      HCC_17 -24.7539498  -2.7295922 HCC_17 HCC_17        17   HCC   Male  57
      HCC_18  -3.5792899   4.2496894 HCC_18 HCC_18        18   HCC   Male  66
      HCC_19 -36.4073373 -18.9035671 HCC_19 HCC_19        19   HCC Female  56
      C_1     17.8693618 -27.2894657    C_1    C_1         1     C Female  57
      C_2     11.0231420  10.3604848    C_2    C_2         2     C   Male  71
      C_3     15.9548510  10.0346953    C_3    C_3         3     C Female  42
      C_4      3.0549555   5.7549929    C_4    C_4         4     C   Male  73
      C_5     23.3247927 -16.1928979    C_5    C_5         5     C   Male  51
      C_6     20.9111502 -29.2877506    C_6    C_6         6     C   Male  21
      C_7     13.9851717   7.2886426    C_7    C_7         7     C   Male  71
      C_8     16.4244088 -20.9927388    C_8    C_8         8     C   Male  74
      C_9     29.1398235 -29.3764869    C_9    C_9         9     C   Male  65
      C_10     9.9247225   9.8194045   C_10   C_10        10     C Female  59
      C_11     9.8373391   9.6411225   C_11   C_11        11     C Female  76
      C_12    20.2163418  10.6614205   C_12   C_12        12     C   Male  53
      C_13    11.5776093   8.3909430   C_13   C_13        13     C   Male  62
      C_14     3.2860558   8.9931948   C_14   C_14        14     C   Male  72
      C_15    12.7079189  10.0911443   C_15   C_15        15     C Female  67
      C_16    18.7842415  -4.3754468   C_16   C_16        16     C   Male  75
      C_17    14.2440095  10.3886285   C_17   C_17        17     C   Male  57
      C_18    12.9278292   8.2691527   C_18   C_18        18     C   Male  66
      C_19     9.9243299   7.5171980   C_19   C_19        19     C Female  56
             Tumor_Staging Lymph_Nodes Grading Bloodvessel_status Resection_Margin
      HCC_1            pT1         pN0      G3                 V0               R0
      HCC_2            pT2         pNX      G2                 V0               R0
      HCC_3            pT1         pNX      G3                 V0               R0
      HCC_4            pT1          NX      G2                 V0               R0
      HCC_5            pT1         pN0      G1                 V0               R0
      HCC_6            pT2        pN1       G2                 V1               R0
      HCC_7            pT1         pNX      G3                 V0               R0
      HCC_8            pT1         pNX      G1                 V0               R0
      HCC_9            pT1          NX      G3                 V0               R0
      HCC_10          pT3a          NX      G2                 V0               R0
      HCC_11          pT2a          NX      G3                 V1               R0
      HCC_12          pT3a          NX      G2                 V1               R0
      HCC_13          pT2a          NX      G2                 V1               R0
      HCC_14           pT1         NA       G1                NA                R1
      HCC_15           pT1          N0      G1                 V0               R0
      HCC_16           pT1          N0      G1                 V0               R0
      HCC_17           pT2         pNX      G2                 V0               R0
      HCC_18           pT1         pNX      G2                 V0               R0
      HCC_19          pT3a         pNX      G3                 V0               R0
      C_1             <NA>        <NA>    <NA>               <NA>             <NA>
      C_2             <NA>        <NA>    <NA>               <NA>             <NA>
      C_3             <NA>        <NA>    <NA>               <NA>             <NA>
      C_4             <NA>        <NA>    <NA>               <NA>             <NA>
      C_5             <NA>        <NA>    <NA>               <NA>             <NA>
      C_6             <NA>        <NA>    <NA>               <NA>             <NA>
      C_7             <NA>        <NA>    <NA>               <NA>             <NA>
      C_8             <NA>        <NA>    <NA>               <NA>             <NA>
      C_9             <NA>        <NA>    <NA>               <NA>             <NA>
      C_10            <NA>        <NA>    <NA>               <NA>             <NA>
      C_11            <NA>        <NA>    <NA>               <NA>             <NA>
      C_12            <NA>        <NA>    <NA>               <NA>             <NA>
      C_13            <NA>        <NA>    <NA>               <NA>             <NA>
      C_14            <NA>        <NA>    <NA>               <NA>             <NA>
      C_15            <NA>        <NA>    <NA>               <NA>             <NA>
      C_16            <NA>        <NA>    <NA>               <NA>             <NA>
      C_17            <NA>        <NA>    <NA>               <NA>             <NA>
      C_18            <NA>        <NA>    <NA>               <NA>             <NA>
      C_19            <NA>        <NA>    <NA>               <NA>             <NA>
             Underlying_Liver_Disease Sample
      HCC_1                      NASH  HCC_1
      HCC_2                      NASH  HCC_2
      HCC_3                      <NA>  HCC_3
      HCC_4                      HepC  HCC_4
      HCC_5                       NA   HCC_5
      HCC_6                        NA  HCC_6
      HCC_7                      NASH  HCC_7
      HCC_8                        NA  HCC_8
      HCC_9                      NASH  HCC_9
      HCC_10                       NA HCC_10
      HCC_11                       NA HCC_11
      HCC_12                     HepC HCC_12
      HCC_13                     HepC HCC_13
      HCC_14                HepB/NASH HCC_14
      HCC_15                     NASH HCC_15
      HCC_16    Hep B/Liver cirrhosis HCC_16
      HCC_17     Liver cirrhosis/NASH HCC_17
      HCC_18          Liver cirrhosis HCC_18
      HCC_19                      NA  HCC_19
      C_1                        <NA>    C_1
      C_2                        <NA>    C_2
      C_3                        <NA>    C_3
      C_4                        <NA>    C_4
      C_5                        <NA>    C_5
      C_6                        <NA>    C_6
      C_7                        <NA>    C_7
      C_8                        <NA>    C_8
      C_9                        <NA>    C_9
      C_10                       <NA>   C_10
      C_11                       <NA>   C_11
      C_12                       <NA>   C_12
      C_13                       <NA>   C_13
      C_14                       <NA>   C_14
      C_15                       <NA>   C_15
      C_16                       <NA>   C_16
      C_17                       <NA>   C_17
      C_18                       <NA>   C_18
      C_19                       <NA>   C_19

---

    Code
      PCA$filtered_data
    Output
      NULL

# Test PCA plot with only shape

    Code
      PCA$D_PCA_plot
    Output
                     PCx         PCy  label Sample PatientID Group Gender Age
      HCC_1  -24.2260058  -5.8796975  HCC_1  HCC_1         1   HCC Female  57
      HCC_2   -0.2701998  11.7565454  HCC_2  HCC_2         2   HCC   Male  71
      HCC_3  -32.8285745   2.9649821  HCC_3  HCC_3         3   HCC Female  42
      HCC_4   -8.0581633   3.2113162  HCC_4  HCC_4         4   HCC   Male  73
      HCC_5    5.4800309   4.1258811  HCC_5  HCC_5         5   HCC   Male  51
      HCC_6   -9.0077163  -0.3124016  HCC_6  HCC_6         6   HCC   Male  21
      HCC_7  -16.6865955  -0.9581829  HCC_7  HCC_7         7   HCC   Male  71
      HCC_8   -0.7120405   7.0650912  HCC_8  HCC_8         8   HCC   Male  74
      HCC_9   -8.8132365   4.0503016  HCC_9  HCC_9         9   HCC   Male  65
      HCC_10 -46.4426474 -16.5206031 HCC_10 HCC_10        10   HCC Female  59
      HCC_11 -13.7105749   2.8033642 HCC_11 HCC_11        11   HCC Female  76
      HCC_12 -14.0413060   0.3658427 HCC_12 HCC_12        12   HCC   Male  53
      HCC_13 -14.5881346   3.4584455 HCC_13 HCC_13        13   HCC   Male  62
      HCC_14   3.3520552   6.3015778 HCC_14 HCC_14        14   HCC   Male  72
      HCC_15 -19.2159911   2.4651543 HCC_15 HCC_15        15   HCC Female  67
      HCC_16 -10.6083777   2.7896151 HCC_16 HCC_16        16   HCC   Male  75
      HCC_17 -24.7539498  -2.7295922 HCC_17 HCC_17        17   HCC   Male  57
      HCC_18  -3.5792899   4.2496894 HCC_18 HCC_18        18   HCC   Male  66
      HCC_19 -36.4073373 -18.9035671 HCC_19 HCC_19        19   HCC Female  56
      C_1     17.8693618 -27.2894657    C_1    C_1         1     C Female  57
      C_2     11.0231420  10.3604848    C_2    C_2         2     C   Male  71
      C_3     15.9548510  10.0346953    C_3    C_3         3     C Female  42
      C_4      3.0549555   5.7549929    C_4    C_4         4     C   Male  73
      C_5     23.3247927 -16.1928979    C_5    C_5         5     C   Male  51
      C_6     20.9111502 -29.2877506    C_6    C_6         6     C   Male  21
      C_7     13.9851717   7.2886426    C_7    C_7         7     C   Male  71
      C_8     16.4244088 -20.9927388    C_8    C_8         8     C   Male  74
      C_9     29.1398235 -29.3764869    C_9    C_9         9     C   Male  65
      C_10     9.9247225   9.8194045   C_10   C_10        10     C Female  59
      C_11     9.8373391   9.6411225   C_11   C_11        11     C Female  76
      C_12    20.2163418  10.6614205   C_12   C_12        12     C   Male  53
      C_13    11.5776093   8.3909430   C_13   C_13        13     C   Male  62
      C_14     3.2860558   8.9931948   C_14   C_14        14     C   Male  72
      C_15    12.7079189  10.0911443   C_15   C_15        15     C Female  67
      C_16    18.7842415  -4.3754468   C_16   C_16        16     C   Male  75
      C_17    14.2440095  10.3886285   C_17   C_17        17     C   Male  57
      C_18    12.9278292   8.2691527   C_18   C_18        18     C   Male  66
      C_19     9.9243299   7.5171980   C_19   C_19        19     C Female  56
             Tumor_Staging Lymph_Nodes Grading Bloodvessel_status Resection_Margin
      HCC_1            pT1         pN0      G3                 V0               R0
      HCC_2            pT2         pNX      G2                 V0               R0
      HCC_3            pT1         pNX      G3                 V0               R0
      HCC_4            pT1          NX      G2                 V0               R0
      HCC_5            pT1         pN0      G1                 V0               R0
      HCC_6            pT2        pN1       G2                 V1               R0
      HCC_7            pT1         pNX      G3                 V0               R0
      HCC_8            pT1         pNX      G1                 V0               R0
      HCC_9            pT1          NX      G3                 V0               R0
      HCC_10          pT3a          NX      G2                 V0               R0
      HCC_11          pT2a          NX      G3                 V1               R0
      HCC_12          pT3a          NX      G2                 V1               R0
      HCC_13          pT2a          NX      G2                 V1               R0
      HCC_14           pT1         NA       G1                NA                R1
      HCC_15           pT1          N0      G1                 V0               R0
      HCC_16           pT1          N0      G1                 V0               R0
      HCC_17           pT2         pNX      G2                 V0               R0
      HCC_18           pT1         pNX      G2                 V0               R0
      HCC_19          pT3a         pNX      G3                 V0               R0
      C_1             <NA>        <NA>    <NA>               <NA>             <NA>
      C_2             <NA>        <NA>    <NA>               <NA>             <NA>
      C_3             <NA>        <NA>    <NA>               <NA>             <NA>
      C_4             <NA>        <NA>    <NA>               <NA>             <NA>
      C_5             <NA>        <NA>    <NA>               <NA>             <NA>
      C_6             <NA>        <NA>    <NA>               <NA>             <NA>
      C_7             <NA>        <NA>    <NA>               <NA>             <NA>
      C_8             <NA>        <NA>    <NA>               <NA>             <NA>
      C_9             <NA>        <NA>    <NA>               <NA>             <NA>
      C_10            <NA>        <NA>    <NA>               <NA>             <NA>
      C_11            <NA>        <NA>    <NA>               <NA>             <NA>
      C_12            <NA>        <NA>    <NA>               <NA>             <NA>
      C_13            <NA>        <NA>    <NA>               <NA>             <NA>
      C_14            <NA>        <NA>    <NA>               <NA>             <NA>
      C_15            <NA>        <NA>    <NA>               <NA>             <NA>
      C_16            <NA>        <NA>    <NA>               <NA>             <NA>
      C_17            <NA>        <NA>    <NA>               <NA>             <NA>
      C_18            <NA>        <NA>    <NA>               <NA>             <NA>
      C_19            <NA>        <NA>    <NA>               <NA>             <NA>
             Underlying_Liver_Disease Sample
      HCC_1                      NASH  HCC_1
      HCC_2                      NASH  HCC_2
      HCC_3                      <NA>  HCC_3
      HCC_4                      HepC  HCC_4
      HCC_5                       NA   HCC_5
      HCC_6                        NA  HCC_6
      HCC_7                      NASH  HCC_7
      HCC_8                        NA  HCC_8
      HCC_9                      NASH  HCC_9
      HCC_10                       NA HCC_10
      HCC_11                       NA HCC_11
      HCC_12                     HepC HCC_12
      HCC_13                     HepC HCC_13
      HCC_14                HepB/NASH HCC_14
      HCC_15                     NASH HCC_15
      HCC_16    Hep B/Liver cirrhosis HCC_16
      HCC_17     Liver cirrhosis/NASH HCC_17
      HCC_18          Liver cirrhosis HCC_18
      HCC_19                      NA  HCC_19
      C_1                        <NA>    C_1
      C_2                        <NA>    C_2
      C_3                        <NA>    C_3
      C_4                        <NA>    C_4
      C_5                        <NA>    C_5
      C_6                        <NA>    C_6
      C_7                        <NA>    C_7
      C_8                        <NA>    C_8
      C_9                        <NA>    C_9
      C_10                       <NA>   C_10
      C_11                       <NA>   C_11
      C_12                       <NA>   C_12
      C_13                       <NA>   C_13
      C_14                       <NA>   C_14
      C_15                       <NA>   C_15
      C_16                       <NA>   C_16
      C_17                       <NA>   C_17
      C_18                       <NA>   C_18
      C_19                       <NA>   C_19

---

    Code
      PCA$filtered_data
    Output
      NULL

# Test PCA plot with only colour and continuous colour

    Code
      PCA$D_PCA_plot
    Output
                     PCx         PCy  label Sample PatientID Group Gender Age
      HCC_1  -24.2260058  -5.8796975  HCC_1  HCC_1         1   HCC Female  57
      HCC_2   -0.2701998  11.7565454  HCC_2  HCC_2         2   HCC   Male  71
      HCC_3  -32.8285745   2.9649821  HCC_3  HCC_3         3   HCC Female  42
      HCC_4   -8.0581633   3.2113162  HCC_4  HCC_4         4   HCC   Male  73
      HCC_5    5.4800309   4.1258811  HCC_5  HCC_5         5   HCC   Male  51
      HCC_6   -9.0077163  -0.3124016  HCC_6  HCC_6         6   HCC   Male  21
      HCC_7  -16.6865955  -0.9581829  HCC_7  HCC_7         7   HCC   Male  71
      HCC_8   -0.7120405   7.0650912  HCC_8  HCC_8         8   HCC   Male  74
      HCC_9   -8.8132365   4.0503016  HCC_9  HCC_9         9   HCC   Male  65
      HCC_10 -46.4426474 -16.5206031 HCC_10 HCC_10        10   HCC Female  59
      HCC_11 -13.7105749   2.8033642 HCC_11 HCC_11        11   HCC Female  76
      HCC_12 -14.0413060   0.3658427 HCC_12 HCC_12        12   HCC   Male  53
      HCC_13 -14.5881346   3.4584455 HCC_13 HCC_13        13   HCC   Male  62
      HCC_14   3.3520552   6.3015778 HCC_14 HCC_14        14   HCC   Male  72
      HCC_15 -19.2159911   2.4651543 HCC_15 HCC_15        15   HCC Female  67
      HCC_16 -10.6083777   2.7896151 HCC_16 HCC_16        16   HCC   Male  75
      HCC_17 -24.7539498  -2.7295922 HCC_17 HCC_17        17   HCC   Male  57
      HCC_18  -3.5792899   4.2496894 HCC_18 HCC_18        18   HCC   Male  66
      HCC_19 -36.4073373 -18.9035671 HCC_19 HCC_19        19   HCC Female  56
      C_1     17.8693618 -27.2894657    C_1    C_1         1     C Female  57
      C_2     11.0231420  10.3604848    C_2    C_2         2     C   Male  71
      C_3     15.9548510  10.0346953    C_3    C_3         3     C Female  42
      C_4      3.0549555   5.7549929    C_4    C_4         4     C   Male  73
      C_5     23.3247927 -16.1928979    C_5    C_5         5     C   Male  51
      C_6     20.9111502 -29.2877506    C_6    C_6         6     C   Male  21
      C_7     13.9851717   7.2886426    C_7    C_7         7     C   Male  71
      C_8     16.4244088 -20.9927388    C_8    C_8         8     C   Male  74
      C_9     29.1398235 -29.3764869    C_9    C_9         9     C   Male  65
      C_10     9.9247225   9.8194045   C_10   C_10        10     C Female  59
      C_11     9.8373391   9.6411225   C_11   C_11        11     C Female  76
      C_12    20.2163418  10.6614205   C_12   C_12        12     C   Male  53
      C_13    11.5776093   8.3909430   C_13   C_13        13     C   Male  62
      C_14     3.2860558   8.9931948   C_14   C_14        14     C   Male  72
      C_15    12.7079189  10.0911443   C_15   C_15        15     C Female  67
      C_16    18.7842415  -4.3754468   C_16   C_16        16     C   Male  75
      C_17    14.2440095  10.3886285   C_17   C_17        17     C   Male  57
      C_18    12.9278292   8.2691527   C_18   C_18        18     C   Male  66
      C_19     9.9243299   7.5171980   C_19   C_19        19     C Female  56
             Tumor_Staging Lymph_Nodes Grading Bloodvessel_status Resection_Margin
      HCC_1            pT1         pN0      G3                 V0               R0
      HCC_2            pT2         pNX      G2                 V0               R0
      HCC_3            pT1         pNX      G3                 V0               R0
      HCC_4            pT1          NX      G2                 V0               R0
      HCC_5            pT1         pN0      G1                 V0               R0
      HCC_6            pT2        pN1       G2                 V1               R0
      HCC_7            pT1         pNX      G3                 V0               R0
      HCC_8            pT1         pNX      G1                 V0               R0
      HCC_9            pT1          NX      G3                 V0               R0
      HCC_10          pT3a          NX      G2                 V0               R0
      HCC_11          pT2a          NX      G3                 V1               R0
      HCC_12          pT3a          NX      G2                 V1               R0
      HCC_13          pT2a          NX      G2                 V1               R0
      HCC_14           pT1         NA       G1                NA                R1
      HCC_15           pT1          N0      G1                 V0               R0
      HCC_16           pT1          N0      G1                 V0               R0
      HCC_17           pT2         pNX      G2                 V0               R0
      HCC_18           pT1         pNX      G2                 V0               R0
      HCC_19          pT3a         pNX      G3                 V0               R0
      C_1             <NA>        <NA>    <NA>               <NA>             <NA>
      C_2             <NA>        <NA>    <NA>               <NA>             <NA>
      C_3             <NA>        <NA>    <NA>               <NA>             <NA>
      C_4             <NA>        <NA>    <NA>               <NA>             <NA>
      C_5             <NA>        <NA>    <NA>               <NA>             <NA>
      C_6             <NA>        <NA>    <NA>               <NA>             <NA>
      C_7             <NA>        <NA>    <NA>               <NA>             <NA>
      C_8             <NA>        <NA>    <NA>               <NA>             <NA>
      C_9             <NA>        <NA>    <NA>               <NA>             <NA>
      C_10            <NA>        <NA>    <NA>               <NA>             <NA>
      C_11            <NA>        <NA>    <NA>               <NA>             <NA>
      C_12            <NA>        <NA>    <NA>               <NA>             <NA>
      C_13            <NA>        <NA>    <NA>               <NA>             <NA>
      C_14            <NA>        <NA>    <NA>               <NA>             <NA>
      C_15            <NA>        <NA>    <NA>               <NA>             <NA>
      C_16            <NA>        <NA>    <NA>               <NA>             <NA>
      C_17            <NA>        <NA>    <NA>               <NA>             <NA>
      C_18            <NA>        <NA>    <NA>               <NA>             <NA>
      C_19            <NA>        <NA>    <NA>               <NA>             <NA>
             Underlying_Liver_Disease Sample
      HCC_1                      NASH  HCC_1
      HCC_2                      NASH  HCC_2
      HCC_3                      <NA>  HCC_3
      HCC_4                      HepC  HCC_4
      HCC_5                       NA   HCC_5
      HCC_6                        NA  HCC_6
      HCC_7                      NASH  HCC_7
      HCC_8                        NA  HCC_8
      HCC_9                      NASH  HCC_9
      HCC_10                       NA HCC_10
      HCC_11                       NA HCC_11
      HCC_12                     HepC HCC_12
      HCC_13                     HepC HCC_13
      HCC_14                HepB/NASH HCC_14
      HCC_15                     NASH HCC_15
      HCC_16    Hep B/Liver cirrhosis HCC_16
      HCC_17     Liver cirrhosis/NASH HCC_17
      HCC_18          Liver cirrhosis HCC_18
      HCC_19                      NA  HCC_19
      C_1                        <NA>    C_1
      C_2                        <NA>    C_2
      C_3                        <NA>    C_3
      C_4                        <NA>    C_4
      C_5                        <NA>    C_5
      C_6                        <NA>    C_6
      C_7                        <NA>    C_7
      C_8                        <NA>    C_8
      C_9                        <NA>    C_9
      C_10                       <NA>   C_10
      C_11                       <NA>   C_11
      C_12                       <NA>   C_12
      C_13                       <NA>   C_13
      C_14                       <NA>   C_14
      C_15                       <NA>   C_15
      C_16                       <NA>   C_16
      C_17                       <NA>   C_17
      C_18                       <NA>   C_18
      C_19                       <NA>   C_19

---

    Code
      PCA$filtered_data
    Output
      NULL

