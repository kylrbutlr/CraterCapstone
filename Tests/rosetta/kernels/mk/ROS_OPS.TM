KPL/MK

Meta-kernel for the Rosetta dataset release v2.1.0
==========================================================================

   This meta-kernel lists the Rosetta SPICE kernels providing coverage
   for the full mission. All of the kernels listed below are stored in
   the Rosetta SPICE Dataset. This set of files and the order in which
   they are listed were picked to provide the best available data and
   the most complete coverage based on the information about the kernels
   available at the time this meta-kernel was produced.


Usage of the Meta-kernel
-------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make use
   of this frame kernel must "load" the kernel normally during program
   initialization. Loading the kernel associates the data items with
   their names in a data structure called the "kernel pool". The SPICELIB
   routine FURNSH loads a kernel into the pool.

   Reading the comments of the binary SPK and C-Kernels is very helpful to
   understand the implemented S/C position and orientation. You can use the
   following NAIF command line utility to extract that information from the
   kernels:

      > commnt -r RATT_*

   The kernels listed below can be obtained from the ESA SPICE FTP server:

      ftp://spiftp.esac.esa.int/data/SPICE/ROSETTA/kernels/


Implementation Notes
-------------------------------------------------------------------------

   It is recommended that users make a local copy of this file and
   modify the value of the PATH_VALUES keyword to point to the actual
   location of the ROSETTA SPICE data set's ``data'' directory on
   their system. Replacing ``/'' with ``\'' and converting line
   terminators to the format native to the user's system may also be
   required if this meta-kernel is to be used on a non-UNIX workstation.

-------------------

   This file was created on October 16, 2018 by Marc Costa Sitja ESAC/ESA.


   \begindata

     PATH_VALUES       = ( '../Tests/rosetta/kernels'  )

     PATH_SYMBOLS      = ( 'KERNELS' )

     KERNELS_TO_LOAD   = (

                          '$KERNELS/lsk/NAIF0011.TLS'


                          '$KERNELS/sclk/ROS_160929_STEP.TSC'

                          '$KERNELS/sclk/LANDER_170904_STEP.TSC'


                          '$KERNELS/pck/PCK00010.TPC'
                          '$KERNELS/pck/DE403_MASSES.TPC'
                          '$KERNELS/pck/former_versions/EARTH_000101_170227_161206.BPC'

                          '$KERNELS/pck/ROS_LUTETIA_RSOC_V03.TPC'
                          '$KERNELS/pck/ROS_STEINS_V05.TPC'
                          '$KERNELS/pck/ROS_CG_RAD_V10.TPC'


                          '$KERNELS/fk/EARTHFIXEDITRF93.TF'
                          '$KERNELS/fk/EARTH_TOPO_050714.TF'
                          '$KERNELS/fk/ESTRACK_V01.TF'
                          '$KERNELS/fk/NEW_NORCIA_TOPO.TF'

                          '$KERNELS/fk/ROS_LUTETIA_RSOC_V03.TF'

                          '$KERNELS/fk/RSSD0002.TF'

                          '$KERNELS/fk/ROS_V33.TF'


                          '$KERNELS/ik/ROS_ALICE_V17.TI'
                          '$KERNELS/ik/ROS_CIVA_V10.TI'
                          '$KERNELS/ik/ROS_CONSERT_V11.TI'
                          '$KERNELS/ik/ROS_COSIMA_V14.TI'
                          '$KERNELS/ik/ROS_DIM_V11.TI'
                          '$KERNELS/ik/ROS_GIADA_V12.TI'
                          '$KERNELS/ik/ROS_MIDAS_V11.TI'
                          '$KERNELS/ik/ROS_MIRO_V11.TI'
                          '$KERNELS/ik/ROS_NAVCAM_V03.TI'
                          '$KERNELS/ik/ROS_OSIRIS_V15.TI'
                          '$KERNELS/ik/ROS_ROLIS_V03.TI'
                          '$KERNELS/ik/ROS_ROSINA_V12.TI'
                          '$KERNELS/ik/ROS_RPC_V19.TI'
                          '$KERNELS/ik/ROS_SREM_V01.TI'
                          '$KERNELS/ik/ROS_STR_V11.TI'
                          '$KERNELS/ik/ROS_VIRTIS_V14.TI'


                          '$KERNELS/spk/DE405.BSP'
                          '$KERNELS/spk/MAR097_030101_170101_V0001.BSP'
                          '$KERNELS/spk/OUTERPLANETS_V0003.BSP'

                          '$KERNELS/spk/EARTHSTNS_ITRF93_050714.BSP'
                          '$KERNELS/spk/ESTRACK_V01.BSP'
                          '$KERNELS/spk/NEW_NORCIA.BSP'

                          '$KERNELS/spk/TEMPEL1_9P_DI_P.BSP'
                          '$KERNELS/spk/VESTA_HOR_V01.BSP'
                          '$KERNELS/spk/C2010_A2_LINEAR_HOR_V01.BSP'
                          '$KERNELS/spk/C2004_Q2_MACHHOLZ_HOR_V01.BSP'
                          '$KERNELS/spk/C2002_T7_LINEAR_HOR_V01.BSP'

                          '$KERNELS/spk/former_versions/21_LUTETIA_2004_2016.BSP'
                          '$KERNELS/spk/ORHS_______________00109.BSP'

                          '$KERNELS/spk/former_versions/2867_STEINS_2004_2016.BSP'
                          '$KERNELS/spk/ORHO_______________00077.BSP'

                          '$KERNELS/spk/former_versions/ORHW_______________00016.BSP'
                          '$KERNELS/spk/ORHW_______________00122.BSP'
                          '$KERNELS/spk/CORB_DV_257_03___T19_00345.BSP'

                          '$KERNELS/spk/ROS_COG_V1.BSP'
                          '$KERNELS/spk/ROS_STRUCT_V6.BSP'

                          '$KERNELS/spk/LORB_ROS_SC_PRESEP_V1_0.BSP'
                          '$KERNELS/spk/LORB_SUN_J2000_SDL_V1_1.BSP'
                          '$KERNELS/spk/SPICE_PHILAE_CFF_SONC_V2_0.BSP'
                          '$KERNELS/spk/LORB_SUN_J2000_RBD_1_V1_2.BSP'

                          '$KERNELS/spk/ORHR___________T19_00122.BSP'
                          '$KERNELS/spk/RORB_DV_257_03___T19_00345.BSP'
                          '$KERNELS/spk/ROS_ORBITER_EXTENSION_V2.BSP'


                          '$KERNELS/ck/ROS_LBOOM_V0.BC'

                          '$KERNELS/ck/ROS_HGA_2004_V0001.BC'
                          '$KERNELS/ck/ROS_HGA_2005_V0001.BC'
                          '$KERNELS/ck/ROS_HGA_2006_V0001.BC'
                          '$KERNELS/ck/ROS_HGA_2007_V0001.BC'
                          '$KERNELS/ck/ROS_HGA_2008_V0019.BC'
                          '$KERNELS/ck/ROS_HGA_2009_V0052.BC'
                          '$KERNELS/ck/ROS_HGA_2010_V0046.BC'
                          '$KERNELS/ck/ROS_HGA_2011_V0018.BC'
                          '$KERNELS/ck/ROS_HGA_2014_V0044.BC'
                          '$KERNELS/ck/ROS_HGA_2015_V0053.BC'
                          '$KERNELS/ck/ROS_HGA_2016_V0042.BC'

                          '$KERNELS/ck/ROS_SA_2004_V0001.BC'
                          '$KERNELS/ck/ROS_SA_2005_V0001.BC'
                          '$KERNELS/ck/ROS_SA_2006_V0001.BC'
                          '$KERNELS/ck/ROS_SA_2007_V0001.BC'
                          '$KERNELS/ck/ROS_SA_2008_V0038.BC'
                          '$KERNELS/ck/ROS_SA_2009_V0054.BC'
                          '$KERNELS/ck/ROS_SA_2010_V0052.BC'
                          '$KERNELS/ck/ROS_SA_2011_V0020.BC'
                          '$KERNELS/ck/ROS_SA_2014_V0047.BC'
                          '$KERNELS/ck/ROS_SA_2015_V0042.BC'
                          '$KERNELS/ck/ROS_SA_2016_V0041.BC'

                          '$KERNELS/ck/ROS_VIRTIS_M_0403_1509_V1.BC'
                          '$KERNELS/ck/ROS_VIRTIS_M_1509_1609_V1.BC'

                          '$KERNELS/ck/ATNR_P040302093352_T6_00127.BC'
                          '$KERNELS/ck/former_versions/RATT_DV_145_01_01_T6_00216.BC'
                          '$KERNELS/ck/former_versions/RATT_DV_223_01_01_T6_00302.BC'
                          '$KERNELS/ck/RATT_DV_257_02_01_T6_00344.BC'
                          '$KERNELS/ck/ROS_ORBITER_EXTENSION_V2.BC'


                          '$KERNELS/ck/former_versions/CATT_DV_145_02_______00216.BC'
                          '$KERNELS/ck/former_versions/CATT_DV_223_02_______00302.BC'
                          '$KERNELS/ck/CATT_DV_257_03_______00344.BC'

                          '$KERNELS/ck/LATT_ROS2LDR_PRESEP_V1_0.BC'
                          '$KERNELS/ck/LATT_EME2LDR_SDL_V1_0.BC'
                          '$KERNELS/ck/LATT_ROS2LDR_SDL_ROMAP_V1_0.BC'
                          '$KERNELS/ck/LATT_CFF2LDR_FSS_V2_0.BC'

                         )

   \begintext


Contact Information
------------------------------------------------------------------------

   If you have any questions regarding this file contact SPICE support at
   ESAC:

           Marc Costa Sitja
           (+34) 91-8131-457
           mcosta@sciops.esa.int, esa_spice@sciops.esa.int


   or NAIF at JPL:

           Boris Semenov
           +1 (818) 354-8136
           Boris.Semenov@jpl.nasa.gov


End of MK file.