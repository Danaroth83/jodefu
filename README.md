# JoDeFu

Implementation in MATLAB of the "multiresolution compressed acquisition" (MRCA) image formation and "joint demosaicing and fusion" (JoDeFu) image reconstruction methods.

The code is able to:
- model an optical acquisition device based on color filter arrays (CFAs) and/or multiresolution sensors;
- estimate an image datacube from their acquisition;
- compare the results of the estimated products with respect to classical demosaicing and sharpening algorithms.

## Demonstrative scripts

The repository contains demo scripts for the experiments provided in the associated article, which, for the user convenience, are all located in the `src\main` folder and the results are saved in the `data\output` folder. Specifically we provide:
- **Image formation**: Scripts which test the quality of image reconstructed starting from acquisitions modeled through a variety of image formation methods:
  - `demo_formation_mrca.m`:  Reconstructed through the JoDeFu v1 algorithm;
  - `demo_formation_classic.m`:  Reconstructed with classic noniterative algorithms;
  - `demo_formation_cassi.m`: Reconstructed from CASSI acquisitions both with the CASSI decoder and JoDeFu v1;
  - `demo_formation_software.m`: Obtained with software compression encoder/decoders.
- **Image reconstruction**: Scripts testing various reconstruction algorithm processing a MRCA acquisition:
  - `demo_reconstruction_jodefu.m`: Reconstructed through JoDeFu v1 and v2.
  - `demo_reconstruction_classic.m`: Reconstructed through cascaded classic demosaicing and sharpening algorithms;
- **Parameters' setting**: Scripts testing the settings of the JoDeFu algorithm:
  - `demo_parameters.m`: By varying the regularization parameter, the metric function norm, the total variation linear operator and the blurring diameter.


## License

This project is licensed under the [MIT] License - see the LICENSE.md file for details.

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)