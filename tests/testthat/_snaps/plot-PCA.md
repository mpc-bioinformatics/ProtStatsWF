# Test PCA plot

    Code
      PCA$D_PCA_plot
    Output
                     PCx         PCy group1 group2  label Sample
      HCC_1  -24.2260058  -5.8796975    HCC Female  HCC_1  HCC_1
      HCC_2   -0.2701998  11.7565454    HCC   Male  HCC_2  HCC_2
      HCC_3  -32.8285745   2.9649821    HCC Female  HCC_3  HCC_3
      HCC_4   -8.0581633   3.2113162    HCC   Male  HCC_4  HCC_4
      HCC_5    5.4800309   4.1258811    HCC   Male  HCC_5  HCC_5
      HCC_6   -9.0077163  -0.3124016    HCC   Male  HCC_6  HCC_6
      HCC_7  -16.6865955  -0.9581829    HCC   Male  HCC_7  HCC_7
      HCC_8   -0.7120405   7.0650912    HCC   Male  HCC_8  HCC_8
      HCC_9   -8.8132365   4.0503016    HCC   Male  HCC_9  HCC_9
      HCC_10 -46.4426474 -16.5206031    HCC Female HCC_10 HCC_10
      HCC_11 -13.7105749   2.8033642    HCC Female HCC_11 HCC_11
      HCC_12 -14.0413060   0.3658427    HCC   Male HCC_12 HCC_12
      HCC_13 -14.5881346   3.4584455    HCC   Male HCC_13 HCC_13
      HCC_14   3.3520552   6.3015778    HCC   Male HCC_14 HCC_14
      HCC_15 -19.2159911   2.4651543    HCC Female HCC_15 HCC_15
      HCC_16 -10.6083777   2.7896151    HCC   Male HCC_16 HCC_16
      HCC_17 -24.7539498  -2.7295922    HCC   Male HCC_17 HCC_17
      HCC_18  -3.5792899   4.2496894    HCC   Male HCC_18 HCC_18
      HCC_19 -36.4073373 -18.9035671    HCC Female HCC_19 HCC_19
      C_1     17.8693618 -27.2894657      C Female    C_1    C_1
      C_2     11.0231420  10.3604848      C   Male    C_2    C_2
      C_3     15.9548510  10.0346953      C Female    C_3    C_3
      C_4      3.0549555   5.7549929      C   Male    C_4    C_4
      C_5     23.3247927 -16.1928979      C   Male    C_5    C_5
      C_6     20.9111502 -29.2877506      C   Male    C_6    C_6
      C_7     13.9851717   7.2886426      C   Male    C_7    C_7
      C_8     16.4244088 -20.9927388      C   Male    C_8    C_8
      C_9     29.1398235 -29.3764869      C   Male    C_9    C_9
      C_10     9.9247225   9.8194045      C Female   C_10   C_10
      C_11     9.8373391   9.6411225      C Female   C_11   C_11
      C_12    20.2163418  10.6614205      C   Male   C_12   C_12
      C_13    11.5776093   8.3909430      C   Male   C_13   C_13
      C_14     3.2860558   8.9931948      C   Male   C_14   C_14
      C_15    12.7079189  10.0911443      C Female   C_15   C_15
      C_16    18.7842415  -4.3754468      C   Male   C_16   C_16
      C_17    14.2440095  10.3886285      C   Male   C_17   C_17
      C_18    12.9278292   8.2691527      C   Male   C_18   C_18
      C_19     9.9243299   7.5171980      C Female   C_19   C_19

---

    Code
      PCA$filtered_data
    Output
      class: SummarizedExperiment 
      dim: 1273 38 
      metadata(0):
      assays(2): intensity_norm intensity
      rownames(1273): P00761 A0A087WW87 ... Q9Y6N5 Q9Y6Y8
      rowData names(5): Protein Gene Protein.Length Organism Description
      colnames(38): HCC_1 HCC_2 ... C_18 C_19
      colData names(11): Sample PatientID ... Resection_Margin
        Underlying_Liver_Disease

