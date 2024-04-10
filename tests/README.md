# tests

Directory containing test cases. 
The following test can be found in the directory: 
- [**One layer, spherically-symmetric**](Test_One_Layer_Spherically_Symmetric.mlx): Compares LOV3D Love numbers against love numbers obtained analytically for a uniform spherically-symmetric bode
- [**Enceladus with lateral variations**](Test_Enceladus_Two_Layers_Lateral_Variations.mlx): 3 layer Enceladus model consisting of a rigid core, ocean and ice-shell with lateral variations. Compares LOV3D Love numbers against love numbers obtained using the perturbation method of [Qin et al.](https://doi.org/10.1093/gji/ggu279) and the FEM model of [Berne et al.](https://doi.org/10.1029/2023GL106656). Reproduces Figure 2 of [Rovira et al. 2024](https://doi.org/10.48550/arXiv.2311.15710)
- [**Multi-layered Spherically-symmetric **](Test_Io_Multi_Layer_Spherically_Symmetric.mlx): Multi-layered Io model based on [Steinke et al. 2020](https://doi.org/10.1016/j.icarus.2019.05.001), consisting of core, deep mantle, asthenosphere and lithosphere. The script obtains the Love numbers and compares them against results obtained with the spherically-symmetric code of  [Rovira-Navarro et al. 2022](https://doi.org/10.1029/2021JE007117)
- [**Multi-layered spherically-symmetric icy moon**](Test_Europa_Titan_Spherically_Symmetric.mlx): Multi-layered icy moon model. The script computes the Love numbers for a multi-layered Europa and Titan models based on [Beuthe et al. 2013](https://www.sciencedirect.com/science/article/pii/S0019103512004745?casa_token=xg0XfpmaHT4AAAAA:Qau6ppdURvhX_Vgm_NiDZVwEtERNnqcosVviHYGaLIHJLBugG7ZgBEnHNPG921Qc5SZAktQ6kw). 




