/*
ADAPTATION DES SORTIES NETCDF ICHTHYOP POUR EVOL-DEB 2013
 */
package util;

import evolution.Simulation;
import evolution.Population;
import evolution.Poissons;
import java.io.File;
import java.io.IOException;
import ucar.nc2.NetcdfFileWriteable;
import ucar.nc2.Dimension;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayInt;
import evolution.Dataset_EVOL;
      import java.util.List;
import ucar.ma2.InvalidRangeException;

//import maltcms.tools.ArrayTools;
import java.util.ArrayList;
import ucar.ma2.DataType;

import ucar.ma2.Array;
/**
 *
 * @author timbrochier
 */
public class OutpoutManager_netcdf {

    static String basename, filename;
    ArrayDouble Var_array;
    private Dimension timeDim, drifter, x_Dim, y_Dim, drifter_data_indiv;
    static private Dimension drifter_TRACE, timeDim_TRACE;
    int nb_drifter_this_year;  
    int compteur_DataIndiv;   
    
    boolean record_OPTIONAL_data = false;
    //+ decommenter ligne 230 et 266 =  Record_OPTIONAL_data

    //   private NCDimFactory dimensionFactory;
    /**
     * Object for creating/writing netCDF files.
     */
    private static NetcdfFileWriteable ncOut, ncOut_Biom2D, ncOut_Data_Indiv, ncOut_TRACE;

// -----------------------------------------------------------------------------
    public void openTrajIndiv_NetCDF() throws IOException {

        System.out.println("Creation fichier netcdf year count " + Simulation.compteur_year);

        String middle_filename = "_TrajIndiv_Y";
        filename = makeFileLocation(middle_filename); //"testWrite.nc";

        ncOut = NetcdfFileWriteable.createNew(filename, true);

        ncOut.addGlobalAttribute("_FillValue", Float.NaN);

        nb_drifter_this_year = Population.Pop.size()*10; //Simulation.Pop_max*10;
        drifter = ncOut.addDimension("drifter", nb_drifter_this_year);//Simulation.Pop_max*10);
//        drifter = ncOut.addUnlimitedDimension("drifter");//Simulation.Pop_max*10);

                timeDim = ncOut.addUnlimitedDimension("time");
//        timeDim = ncOut.addDimension("time", Simulation.nb_jours_par_ans);


        Dimension[] dim2 = new Dimension[2];
        dim2[0] = timeDim;
        dim2[1] = drifter;

        ncOut.addVariable("longitude", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("longitude", "units", "degree east");
        ncOut.addVariableAttribute("longitude", "_FillValue", Double.NaN);

        ncOut.addVariable("latitude", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("latitude", "units", "degree north");
        ncOut.addVariableAttribute("latitude", "_FillValue", Double.NaN);

        ncOut.addVariable("fish_identification", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("fish_identification", "units", "id_number");
        ncOut.addVariableAttribute("fish_identification", "_FillValue", Double.NaN);

        ncOut.addVariable("fish_age", DataType.INT, dim2);
        ncOut.addVariableAttribute("fish_age", "units", "days");
        ncOut.addVariableAttribute("fish_age", "_FillValue", Double.NaN);

        ncOut.addVariable("Effectif_super_indiv", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("Effectif_super_indiv", "units", "N");
        ncOut.addVariableAttribute("Effectif_super_indiv", "_FillValue", Double.NaN);

        ncOut.addVariable("E", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("E", "units", "joules reserve");
        ncOut.addVariableAttribute("E", "_FillValue", Double.NaN);

        ncOut.addVariable("V", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("V", "units", "volume structurel");
        ncOut.addVariableAttribute("V", "_FillValue", Double.NaN);

        ncOut.addVariable("E_R", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("E_R", "units", "joules reproduction buffer");
        ncOut.addVariableAttribute("E_R", "_FillValue", Double.NaN);

        ncOut.addVariable("temperature", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("temperature", "units", "Celsius");
        ncOut.addVariableAttribute("temperature", "_FillValue", Double.NaN);

        ncOut.addVariable("f", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("f", "units", "food functional response");
        ncOut.addVariableAttribute("f", "_FillValue", Double.NaN);
        
        ncOut.addVariable("M_DD", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("M_DD", "units", "Index de Mortalite de competition pour kinesis");
        ncOut.addVariableAttribute("M_DD", "_FillValue", Double.NaN);        
        
        ncOut.addVariable("Nb_eggs_spawned_dt", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("Nb_eggs_spawned_dt", "units", "Number of eggs spawned each day");
        ncOut.addVariableAttribute("Nb_eggs_spawned_dt", "_FillValue", Double.NaN);

        ncOut.addVariable("E_H", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("E_H", "units", "joules maturation buffer");
        ncOut.addVariableAttribute("E_H", "_FillValue", Double.NaN);

        ncOut.addVariable("dx_kin", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("dx_kin", "units", "x-kinesis movement (m/s)");
        ncOut.addVariableAttribute("dx_kin", "_FillValue", Double.NaN);

        ncOut.addVariable("dy_kin", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("dy_kin", "units", "y-kinesis movement (m/s)");
        ncOut.addVariableAttribute("dy_kin", "_FillValue", Double.NaN);

        ncOut.addVariable("dx_advec", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("dx_advec", "units", "x-advection movement (m/s)");
        ncOut.addVariableAttribute("dx_advec", "_FillValue", Double.NaN);

        ncOut.addVariable("dy_advec", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("dy_advec", "units", "y-advection movement (m/s)");
        ncOut.addVariableAttribute("dy_advec", "_FillValue", Double.NaN);

        
        if (record_OPTIONAL_data){

        ncOut.addVariable("depth", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("depth", "units", "meters");
        ncOut.addVariableAttribute("depth", "_FillValue", Double.NaN);

        ncOut.addVariable("Nb_SI_spawned_dt", DataType.INT, dim2);
        ncOut.addVariableAttribute("Nb_SI_spawned_dt", "units", "Number of SI spawned each day");
        ncOut.addVariableAttribute("Nb_SI_spawned_dt", "_FillValue", Double.NaN);

        ncOut.addVariable("Q", DataType.DOUBLE, dim2);
        ncOut.addVariableAttribute("Q", "units", "Habitats quality for kinesis movement");
        ncOut.addVariableAttribute("Q", "_FillValue", Double.NaN);

        
        }

        // create the file
        try {
            ncOut.create();
        } catch (IOException e) {
            System.err.println("ERROR creating file " + ncOut.getLocation() + "\n" + e);
        }
    }

    public void write_TrajIndiv_NetCDF(ArrayList pop) throws Exception {

            //System.out.println("write_TrajIndiv_NetCDF TrajIndiv...  pop.size = " + pop.size());

//       ncOut = NetcdfFileWriteable.openExisting(filename, true);
        //// heres where we write the record variables

        // different ways to create the data arrays.
        // Note the outer dimension has shape 1, since we will write one record at a timeDim
        ArrayDouble.D2 longitude_Data = new ArrayDouble.D2(1,drifter.getLength());
        ArrayDouble.D2 latitude_Data = new ArrayDouble.D2(1,drifter.getLength());
        ArrayDouble.D2 identification_Data = new ArrayDouble.D2(1,drifter.getLength());
        ArrayInt.D2 age_Data = new ArrayInt.D2(1, drifter.getLength()); // <<<<<<<peut-être viré
        ArrayDouble.D2 Effectif_super_indiv_Data = new ArrayDouble.D2(1, drifter.getLength()); // <<<<<<<peut-être viré

// DEB
        ArrayDouble.D2 E_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 V_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 E_R_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 E_H_Data = new ArrayDouble.D2(1, drifter.getLength());

        // ENVIRONNEMENT (NON ESSENTIEL)
        ArrayDouble.D2 temperature_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 f_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 M_DD_Data = new ArrayDouble.D2(1, drifter.getLength());


        
// POUR TESTER LA PONTE :
        ArrayDouble.D2 Nb_eggs_spawned_dt_Data = new ArrayDouble.D2(1, drifter.getLength()); // <<<<<<<peut-être viré

     //   if (record_OPTIONAL_data){
        ArrayInt.D2 prof_Data = new ArrayInt.D2(1, drifter.getLength());
// POUR TESTER LA KINESIS :
        ArrayDouble.D2 Q_Data = new ArrayDouble.D2(1, drifter.getLength());
        
// Vitesse Nage et Advection
        ArrayDouble.D2 dx_kin_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 dy_kin_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 dx_advec_Data = new ArrayDouble.D2(1, drifter.getLength());
        ArrayDouble.D2 dy_advec_Data = new ArrayDouble.D2(1, drifter.getLength());

        ArrayDouble.D2 Nb_SI_spawned_dt_Data = new ArrayDouble.D2(1, drifter.getLength()); // <<<<<<<peut-être viré

                for (int t = 0; t < drifter.getLength(); t++) {
        Q_Data.set(0,t, Double.NaN);
        dx_kin_Data.set(0,t, Double.NaN);
        dy_kin_Data.set(0,t, Double.NaN);
        dx_advec_Data.set(0,t, Double.NaN);
        dy_advec_Data.set(0,t, Double.NaN);        
        
        Nb_SI_spawned_dt_Data.set(0,t, Double.NaN);
        E_H_Data.set(0,t, Double.NaN);
        }

    //    }


       // FILL VALUES
                for (int t = 0; t < drifter.getLength(); t++) {
                        longitude_Data.set(0,t, Double.NaN);
                        latitude_Data.set(0,t, Double.NaN);
                        age_Data.set(0,t, -999);
                        Effectif_super_indiv_Data.set(0,t, -999);
                        E_Data.set(0,t, Double.NaN);
                        V_Data.set(0,t, Double.NaN);
                        E_R_Data.set(0,t, Double.NaN);
                        temperature_Data.set(0,t, Double.NaN);
                        f_Data.set(0,t, Double.NaN);
                        identification_Data.set(0,t, -999);
                        Nb_eggs_spawned_dt_Data.set(0,t, Double.NaN);
                        M_DD_Data.set(0,t, Double.NaN);                        
                }

//int nb_fish_record = Math.min(pop.size(), nb_drifter_this_year);
            int p=0;
        for (int pp = 0; pp < pop.size(); pp++) {
            
            Poissons Po = (Poissons) pop.get(pp);

//            System.out.println("write_TrajIndiv_NetCDF POISSON " + p + " Po.age =  " + Po.age + " Po.id = " + Po.id);
//(Po.age > 30) && (Po.Body_length > 5)
            if ((Simulation.TEST_DEB_flag) || Po.isRecruited){
                //System.out.println("write_TrajIndiv_NetCDF POISSON " + p + " Po.age =  " + Po.age + " Po.id = " + Po.id);
                // On fabrique les data.set :
                longitude_Data.set(0,p, Po.lon);
                latitude_Data.set(0,p, Po.lat);
                identification_Data.set(0,p, Po.id);

                age_Data.set(0, p, Po.age);
                Effectif_super_indiv_Data.set(0, p, Po.S);
                
                E_Data.set(0, p, Po.DEB.E);
                V_Data.set(0, p, Po.DEB.V);
                E_R_Data.set(0, p, Po.DEB.E_R);
                E_H_Data.set(0, p, Po.DEB.E_H);

                temperature_Data.set(0, p, Po.temp);
                f_Data.set(0, p, Po.DEB.f);
                M_DD_Data.set(0, p, Po.M_DDkinesis);

                Nb_eggs_spawned_dt_Data.set(0, p, Po.Nb_eggs_spawned_dt);

//----- Record_OPTIONAL_data *--------------------------------------------******
                Nb_SI_spawned_dt_Data.set(0, p, Po.Nb_SI_spawned_dt );
                prof_Data.set(0, p, (int) Po.depth);
                // Qualite habitat :
                Q_Data.set(0, p, Po.Q);
                // Vitesse de kinesis
                dx_kin_Data.set(0, p, Po.dx_kinesis_swim_km);     
                dy_kin_Data.set(0, p, Po.dy_kinesis_swim_km);     

                dx_advec_Data.set(0, p, Po.dx_advec_km);     
                dy_advec_Data.set(0, p, Po.dy_advec_km);     
                
                //*
//                 */

                // On donne l'origine : le pas de temps actuel :
                int[] origin = new int[]{0, 0};
                origin[0] = Simulation.compteur_t_rec_intra_year;

                // On ecris les data.set dans le netcdf :
                try {
                    // origin = 0
                    ncOut.write("longitude", origin, longitude_Data);
                    ncOut.write("latitude", origin, latitude_Data);
                    ncOut.write("fish_identification", origin, identification_Data);
                    ncOut.write("fish_age", origin, age_Data);
                    ncOut.write("Effectif_super_indiv", origin, Effectif_super_indiv_Data);

                    ncOut.write("E", origin, E_Data);
                    ncOut.write("V", origin, V_Data);
                    ncOut.write("E_R", origin, E_R_Data);
                    ncOut.write("E_H", origin, E_H_Data);

                    ncOut.write("temperature", origin, temperature_Data);
                    ncOut.write("f", origin, f_Data);
                    ncOut.write("M_DD", origin, M_DD_Data);

                    ncOut.write("Nb_eggs_spawned_dt", origin, Nb_eggs_spawned_dt_Data);

//----- Record_OPTIONAL_data *--------------------------------------------******

                   // ncOut.write("depth", origin, prof_Data);
                //    ncOut.write("Q", origin, Q_Data);
                  //  ncOut.write("Nb_SI_spawned_dt", origin, Nb_SI_spawned_dt_Data);
                    ncOut.write("dx_kin", origin, dx_kin_Data);
                    ncOut.write("dy_kin", origin, dy_kin_Data);
                    ncOut.write("dx_advec", origin, dx_advec_Data);
                    ncOut.write("dy_advec", origin, dy_advec_Data);

                } catch (IOException e) {
                    e.printStackTrace();
                }
                p++;
            }// si (Po.age > 0)
        }// boucle sur pop
    }

    public void CloseResults_NetCDF() throws Exception {

        try {
            ncOut.close();
            String strFilePart = ncOut.getLocation();
            String strFileBase = strFilePart.substring(0, strFilePart.indexOf(".part"));
            File filePart = new File(strFilePart);
            File fileBase = new File(strFileBase);
            filePart.renameTo(fileBase);
            System.out.println("Closed NetCDF output file.");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

// ------------------ Biom2D -------------------------------------
    public void open_Biom2D_NetCDF() throws IOException {

        System.out.println("Creation fichier netcdf BIOM 2D year count " + Simulation.compteur_year);

        String middle_filename = "_Biom2D_Y";

        filename = makeFileLocation(middle_filename); //"testWrite.nc";

        ncOut_Biom2D = NetcdfFileWriteable.createNew(filename, true);
        ncOut_Biom2D.addGlobalAttribute("_FillValue", Float.NaN);

        // add dimensions, including unlimited
        x_Dim = ncOut_Biom2D.addDimension("xi_rho", Dataset_EVOL.get_nx());
        y_Dim = ncOut_Biom2D.addDimension("eta_rho", Dataset_EVOL.get_ny());
        timeDim = ncOut_Biom2D.addUnlimitedDimension("time");

        Dimension[] dim3 = new Dimension[3];
        dim3[0] = timeDim;
        dim3[1] = x_Dim;
        dim3[2] = y_Dim;

/*        ncOut_Biom2D.addVariable("x", DataType.FLOAT, new Dimension[] {x_Dim});
        ncOut_Biom2D.addVariableAttribute("x", "standard_name", "x_grid_index");

        ncOut_Biom2D.addVariable("y", DataType.FLOAT, new Dimension[] {y_Dim});
        ncOut_Biom2D.addVariableAttribute("y", "standard_name", "y_grid_index");
*/
        
    // CREATION DES VARIABLES QUI CONTIENDRONT LES VALEURS DES DIMMENSIONS (lon, lat et time)  

    ncOut_Biom2D.addVariable("xi_rho", DataType.FLOAT, new Dimension[] {x_Dim});
    ncOut_Biom2D.addVariableAttribute("xi_rho",  "units","seconds since 1977-01-01");
    ncOut_Biom2D.addVariableAttribute("xi_rho",  "long_name" ,"y-dimension of the grid");
    ncOut_Biom2D.addVariableAttribute("xi_rho",  "field","time, scalar, series");
    ncOut_Biom2D.addVariableAttribute("xi_rho",  "standard_name" ,"y_grid_index");
    ncOut_Biom2D.addVariableAttribute("xi_rho",  "c_grid_dynamic_range" ,"2:217");
    ncOut_Biom2D.addVariableAttribute("xi_rho",  "axis" ,"X");
            
    ncOut_Biom2D.addVariable("eta_rho", DataType.FLOAT, new Dimension[] {y_Dim});
    ncOut_Biom2D.addVariableAttribute("eta_rho",  "units","seconds since 1977-01-01");
    ncOut_Biom2D.addVariableAttribute("eta_rho",  "long_name" ,"y-dimension of the grid");
    ncOut_Biom2D.addVariableAttribute("eta_rho",  "field","time, scalar, series");
    ncOut_Biom2D.addVariableAttribute("eta_rho",  "standard_name" ,"y_grid_index");
    ncOut_Biom2D.addVariableAttribute("eta_rho",  "c_grid_dynamic_range" ,"2:325");
    ncOut_Biom2D.addVariableAttribute("eta_rho",  "axis" ,"Y");

    ncOut_Biom2D.addVariable("time", DataType.FLOAT, new Dimension[] {timeDim});
    ncOut_Biom2D.addVariableAttribute("time",  "units","seconds since 1977-01-01");
    ncOut_Biom2D.addVariableAttribute("time",  "long_name" ,"averaged time since initialization");
    ncOut_Biom2D.addVariableAttribute("time",  "field","time, scalar, series");
    ncOut_Biom2D.addVariableAttribute("time",  "standard_name" ,"time");
    ncOut_Biom2D.addVariableAttribute("time",  "_CoordinateAxisType" ,"Time");
    ncOut_Biom2D.addVariableAttribute("time",  "axis" ,"T");
//    ncOut_Biom2D.addVariableAttribute("time",  "calendar" ,"julian");


    ncOut_Biom2D.addVariable("lon_rho", DataType.FLOAT, new Dimension[] {x_Dim});
    ncOut_Biom2D.addVariableAttribute("lon_rho",  "units","degrees_east");
    ncOut_Biom2D.addVariableAttribute("lon_rho",  "long_name" ,"longitude of RHO-points");
    ncOut_Biom2D.addVariableAttribute("lon_rho",  "field","lon_rho, scalar, series");
    ncOut_Biom2D.addVariableAttribute("lon_rho",  "standard_name" ,"longitude");
    ncOut_Biom2D.addVariableAttribute("lon_rho",  "axis" ,"X");
    ncOut_Biom2D.addVariableAttribute("lon_rho",  "_CoordinateAxisType" ,"Lon");

    ncOut_Biom2D.addVariable("lat_rho", DataType.FLOAT, new Dimension[] {y_Dim});
    ncOut_Biom2D.addVariableAttribute("lat_rho",  "units","degrees_north");
    ncOut_Biom2D.addVariableAttribute("lat_rho",  "long_name" ,"latitude of RHO-points");
    ncOut_Biom2D.addVariableAttribute("lat_rho",  "field","lat_rho, scalar, series");
    ncOut_Biom2D.addVariableAttribute("lat_rho",  "standard_name" ,"latitude");
    ncOut_Biom2D.addVariableAttribute("lat_rho",  "axis" ,"Y");
    ncOut_Biom2D.addVariableAttribute("lat_rho",  "_CoordinateAxisType" ,"Lat");

    
    // "Vraie" Variables :
        ncOut_Biom2D.addVariable("Biomass_ichthyo", DataType.FLOAT, dim3);
        ncOut_Biom2D.addVariableAttribute("Biomass_ichthyo", "units", "tons");
        ncOut_Biom2D.addVariableAttribute("Biomass_ichthyo", "_FillValue", Float.NaN);

        ncOut_Biom2D.addVariable("Biomass_juvenile", DataType.FLOAT, dim3);
        ncOut_Biom2D.addVariableAttribute("Biomass_juvenile", "units", "tons");
        ncOut_Biom2D.addVariableAttribute("Biomass_juvenile", "_FillValue", Float.NaN);

        ncOut_Biom2D.addVariable("Biomass_moy", DataType.FLOAT, dim3);
        ncOut_Biom2D.addVariableAttribute("Biomass_moy", "units", "tons");
        ncOut_Biom2D.addVariableAttribute("Biomass_moy", "_FillValue", Float.NaN);

        ncOut_Biom2D.addVariable("Biomass_grand", DataType.FLOAT, dim3);
        ncOut_Biom2D.addVariableAttribute("Biomass_grand", "units", "tons");
        ncOut_Biom2D.addVariableAttribute("Biomass_grand", "_FillValue", Float.NaN);
        
        ncOut_Biom2D.addVariable("Biomass_landings", DataType.FLOAT, dim3);
        ncOut_Biom2D.addVariableAttribute("Biomass_landings", "units", "tons");
        ncOut_Biom2D.addVariableAttribute("Biomass_landings", "_FillValue", Float.NaN);
        // create the file
        try {
            ncOut_Biom2D.create();
        } catch (IOException e) {
            System.err.println("ERROR creating file " + ncOut_Biom2D.getLocation() + "\n" + e);
        }
        
// write out the non-record variables
     float[] xi_rho_list = new float[Dataset_EVOL.get_nx()];
     float[] eta_rho_list = new float[Dataset_EVOL.get_ny()];
     
     float[] lonRho_list = new float[Dataset_EVOL.get_nx()];
     float[] latRho_list = new float[Dataset_EVOL.get_ny()];
     
     
        int ii = 2;

     for(int i=0; i<Dataset_EVOL.get_nx(); i++)
    {
        xi_rho_list[i] = ii;
        lonRho_list[i] = (float) Dataset_EVOL.lonRho[1][i]; 
        ii = ii+1;
     }
       ii = 2;
     for(int i=0; i<Dataset_EVOL.get_ny(); i++)
    {
        eta_rho_list[i] = ii;
        latRho_list[i] = (float) Dataset_EVOL.latRho[i][1]; 
        ii = ii+1;
     }

               try {
   ncOut_Biom2D.write("xi_rho", Array.factory(xi_rho_list ));
   ncOut_Biom2D.write("eta_rho", Array.factory(eta_rho_list));        
   ncOut_Biom2D.write("lon_rho", Array.factory(lonRho_list));
   ncOut_Biom2D.write("lat_rho", Array.factory(latRho_list));              
                } catch (InvalidRangeException e) {
                    e.printStackTrace();     
    }
    }

    public void write_Biom2D_NetCDF() throws Exception {

        // different ways to create the data arrays.
        // Note the outer dimension has shape 1, since we will write one record at a timeDim
        ArrayInt.D3 Biomass_ichthyo_Data = new ArrayInt.D3(1,x_Dim.getLength(),y_Dim.getLength());
        ArrayInt.D3 Biomass_juvenile_Data = new ArrayInt.D3(1,x_Dim.getLength(),y_Dim.getLength());
        ArrayInt.D3 Biomass_moy_Data = new ArrayInt.D3(1,x_Dim.getLength(),y_Dim.getLength());
        ArrayInt.D3 Biomass_grand_Data = new ArrayInt.D3(1,x_Dim.getLength(),y_Dim.getLength());
        ArrayInt.D3 Biomass_landings_Data = new ArrayInt.D3(1,x_Dim.getLength(),y_Dim.getLength());

       // On copie les donnee de biomasse depuis Population (cette etape pourrais etre optimisee??)
        // On fabrique les data.set :

            for (int i = 0; i < x_Dim.getLength(); i++) {
                for (int j = 0; j < y_Dim.getLength(); j++) {
                    Biomass_ichthyo_Data.set(0,i,j,Population.Biom_ichthyop_2D[i][j]);
                    Biomass_juvenile_Data.set(0,i,j,Population.Biom_petit_2D[i][j]);
                    Biomass_moy_Data.set(0,i,j,Population.Biom_moy_2D[i][j]);
                    Biomass_grand_Data.set(0,i,j,Population.Biom_grd_2D[i][j]);
                    Biomass_landings_Data.set(0,i,j,Population.Biom_landings_2D[i][j]);
                    }
                }


                // On donne l'origine : le pas de temps actuel :

                int[] origin = new int[]{0, 0, 0};
                origin[0] = Simulation.compteur_t_rec_intra_year;

                // On ecris les data.set dans le netcdf :
                try {
                    // origin = 0
                    ncOut_Biom2D.write("Biomass_ichthyo", origin, Biomass_ichthyo_Data);
                    ncOut_Biom2D.write("Biomass_juvenile", origin, Biomass_juvenile_Data);
                    ncOut_Biom2D.write("Biomass_moy", origin, Biomass_moy_Data);
                    ncOut_Biom2D.write("Biomass_grand", origin, Biomass_grand_Data);
                    ncOut_Biom2D.write("Biomass_landings", origin, Biomass_landings_Data);
                    
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

    public void Close_Biom2D_NetCDF() throws Exception {

        try {
            ncOut_Biom2D.close();
            String strFilePart = ncOut_Biom2D.getLocation();
            String strFileBase = strFilePart.substring(0, strFilePart.indexOf(".part"));
            File filePart = new File(strFilePart);
            File fileBase = new File(strFileBase);
            filePart.renameTo(fileBase);
            System.out.println("Closed NetCDF Biom2D output file.");
        } catch (IOException e) {
            e.printStackTrace();
        }
        //GZipFile()
    }


// -------------------- DataIndiv ------------    
   
    public void open_DataIndiv_NetCDF()
        throws IOException
    {
        System.out.println((new StringBuilder()).append("Creation fichier DATA_INDIV netcdf year count ").append(Simulation.compteur_year).toString());
        String middle_filename = "_DataIndiv_Y";
        filename = makeFileLocation(middle_filename);
        ncOut_Data_Indiv = NetcdfFileWriteable.createNew(filename, true);
        ncOut_Data_Indiv.addGlobalAttribute("_FillValue", Float.valueOf((0.0F / 0.0F)));
        drifter_data_indiv = ncOut_Data_Indiv.addUnlimitedDimension("drifter_data_indiv");
        compteur_DataIndiv = 0;
        Dimension dim1[] = new Dimension[1];
        dim1[0] = drifter_data_indiv;
        ncOut_Data_Indiv.addVariable("fish_identification", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("fish_identification", "units", "id_number");
        ncOut_Data_Indiv.addVariableAttribute("fish_identification", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("longitude_natal", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("longitude_natal", "units", "degree east");
        ncOut_Data_Indiv.addVariableAttribute("longitude_natal", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("latitude_natal", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("latitude_natal", "units", "degree north");
        ncOut_Data_Indiv.addVariableAttribute("latitude_natal", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("jour_natal", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("jour_natal", "units", "jour de l'annee");
        ncOut_Data_Indiv.addVariableAttribute("jour_natal", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("year_natal", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("year_natal", "units", "annee de naissance");
        ncOut_Data_Indiv.addVariableAttribute("year_natal", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("Effectif_SI_natal", DataType.INT, dim1);
        ncOut_Data_Indiv.addVariableAttribute("Effectif_SI_natal", "units", "N");
        ncOut_Data_Indiv.addVariableAttribute("Effectif_SI_natal", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("temperature_natal", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("temperature_natal", "units", "Celsius");
        ncOut_Data_Indiv.addVariableAttribute("temperature_natal", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("longitude_rec", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("longitude_rec", "units", "degree east");
        ncOut_Data_Indiv.addVariableAttribute("longitude_rec", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("latitude_rec", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("latitude_rec", "units", "degree north");
        ncOut_Data_Indiv.addVariableAttribute("latitude_rec", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        ncOut_Data_Indiv.addVariable("Body_length_rec", DataType.INT, dim1);
        ncOut_Data_Indiv.addVariableAttribute("Body_length_rec", "units", "cm");
        ncOut_Data_Indiv.addVariableAttribute("Body_length_rec", "_FillValue", Double.valueOf((0.0D / 0.0D)));

// Tolerances hereditaires : 
        ncOut_Data_Indiv.addVariable("rayon_exploration_temporel", DataType.INT, dim1);
        ncOut_Data_Indiv.addVariableAttribute("rayon_exploration_temporel", "units", "days");
        ncOut_Data_Indiv.addVariableAttribute("rayon_exploration_temporel", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        
        ncOut_Data_Indiv.addVariable("rayon_exploration_spatial", DataType.INT, dim1);
        ncOut_Data_Indiv.addVariableAttribute("rayon_exploration_spatial", "units", "km");
        ncOut_Data_Indiv.addVariableAttribute("rayon_exploration_spatial", "_FillValue", Double.valueOf((0.0D / 0.0D)));

        ncOut_Data_Indiv.addVariable("tolerance_temperature", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("tolerance_temperature", "units", "degCelsius");
        ncOut_Data_Indiv.addVariableAttribute("tolerance_temperature", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        
        ncOut_Data_Indiv.addVariable("tolerance_salinite", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("tolerance_salinite", "units", "days");
        ncOut_Data_Indiv.addVariableAttribute("tolerance_salinite", "_FillValue", Double.valueOf((0.0D / 0.0D)));
        
        ncOut_Data_Indiv.addVariable("sigma_t_swim", DataType.DOUBLE, dim1);
        ncOut_Data_Indiv.addVariableAttribute("sigma_t_swim", "units", "days");
        ncOut_Data_Indiv.addVariableAttribute("sigma_t_swim", "_FillValue", Double.valueOf((0.0D / 0.0D)));

        try
        {
            ncOut_Data_Indiv.create();
        }
        catch(IOException e)
        {
            System.err.println((new StringBuilder()).append("ERROR creating file ").append(ncOut.getLocation()).append("\n").append(e).toString());
        }
    }

    public void write_DataIndiv_NetCDF(Poissons Po)
        throws Exception
    {
        ucar.ma2.ArrayDouble.D1 longitude_natal_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 latitude_natal_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 identification_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 jour_natal_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 year_natal_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 temperature_natal_Data = new ucar.ma2.ArrayDouble.D1(1);

        ucar.ma2.ArrayDouble.D1 rayon_exploration_temporel_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 rayon_exploration_spatial_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 tolerance_temperature_Data = new ucar.ma2.ArrayDouble.D1(1);
        ucar.ma2.ArrayDouble.D1 tolerance_salinite_Data = new ucar.ma2.ArrayDouble.D1(1);
        
        ucar.ma2.ArrayDouble.D1 sigma_t_swim_Data = new ucar.ma2.ArrayDouble.D1(1);
                
        int p = 0;
        identification_Data.set(p, Po.id);
        longitude_natal_Data.set(p, Po.lon_init);
        latitude_natal_Data.set(p, Po.lat_init);
        temperature_natal_Data.set(p, Po.temperature_natal);
        year_natal_Data.set(p, Po.year_init);
        jour_natal_Data.set(p, Po.day_init);

        rayon_exploration_temporel_Data.set(p, Po.rayon_exploration_temporel);
        rayon_exploration_spatial_Data.set(p, Po.rayon_exploration_spatial);
        tolerance_temperature_Data.set(p, Po.tolerance_temperature);
        tolerance_salinite_Data.set(p, Po.tolerance_salinite);
        sigma_t_swim_Data.set(p, Po.sigma_t_swim);
        
        
        int origin[] = {
            0
        };
        origin[0] = compteur_DataIndiv;
        try
        {
            ncOut_Data_Indiv.write("fish_identification", origin, identification_Data);
            ncOut_Data_Indiv.write("longitude_natal", origin, longitude_natal_Data);
            ncOut_Data_Indiv.write("latitude_natal", origin, latitude_natal_Data);
            ncOut_Data_Indiv.write("jour_natal", origin, jour_natal_Data);
            ncOut_Data_Indiv.write("year_natal", origin, year_natal_Data);
            ncOut_Data_Indiv.write("temperature_natal", origin, temperature_natal_Data);

            ncOut_Data_Indiv.write("rayon_exploration_temporel", origin, rayon_exploration_temporel_Data);
            ncOut_Data_Indiv.write("rayon_exploration_spatial", origin, rayon_exploration_spatial_Data);
            ncOut_Data_Indiv.write("tolerance_temperature", origin, tolerance_temperature_Data);
            ncOut_Data_Indiv.write("tolerance_salinite", origin, tolerance_salinite_Data);

            ncOut_Data_Indiv.write("sigma_t_swim", origin, sigma_t_swim_Data);
        }
        catch(IOException e)
        {
            e.printStackTrace();
        }
        compteur_DataIndiv++;
    }

    public void Close_DataIndiv_NetCDF()
        throws Exception
    {
        try
        {
            ncOut_Data_Indiv.close();
            String strFilePart = ncOut_Data_Indiv.getLocation();
            String strFileBase = strFilePart.substring(0, strFilePart.indexOf(".part"));
            File filePart = new File(strFilePart);
            File fileBase = new File(strFileBase);
            filePart.renameTo(fileBase);
            System.out.println("Closed NetCDF output file.");
        }
        catch(IOException e)
        {
            e.printStackTrace();
        }
    }
 
    
    // all done
    // ------------- UTIL ----------------
    private static String makeFileLocation(String middle_filename) throws IOException {

        filename = Simulation.output_dir + Simulation.Strategie + middle_filename + Simulation.compteur_year + ".nc";
        System.out.println(" Création du fichier de sortie NetCDF " + filename);

        File file = new File(filename);
        try {
            //IOTools.makeDirectories(file.getAbsolutePath());
            file.createNewFile();
            file.delete();
        } catch (Exception ex) {
            IOException ioex = new IOException("{Ouput} Failed to create NetCDF file " + filename + " ==> " + ex.getMessage());
            ioex.setStackTrace(ex.getStackTrace());
            throw ioex;
        }
        basename = filename;
        return filename + ".part";
    }


// POUR SAUVEGARDER UNE TRACE (TRAJECTOIRE A HAUTE FREQUENCE, CHAQUE PAS DE TEMPS)
    // ATTENTION, CA PREND DE LA PLACE!!!
// -----------------------------------------------------------------------------
    public static void openTRACE_NetCDF(long fish_id) throws IOException {

        System.out.println("Creation fichier netcdf TRACE - fish_id " + fish_id);

        String middle_filename = "_TRACE_Indiv_ID" + fish_id;
        filename = makeFileLocation(middle_filename); //"testWrite.nc";

        ncOut_TRACE = NetcdfFileWriteable.createNew(filename, true);

        ncOut_TRACE.addGlobalAttribute("_FillValue", Float.NaN);

        drifter_TRACE = ncOut_TRACE.addDimension("drifter", 1);//Simulation.Pop_max*10);

                timeDim_TRACE = ncOut_TRACE.addUnlimitedDimension("time");


        Dimension[] dim2 = new Dimension[2];
        dim2[0] = timeDim_TRACE;
        dim2[1] = drifter_TRACE;

        ncOut_TRACE.addVariable("longitude", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("longitude", "units", "degree east");
        ncOut_TRACE.addVariableAttribute("longitude", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("latitude", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("latitude", "units", "degree north");
        ncOut_TRACE.addVariableAttribute("latitude", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("fish_identification", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("fish_identification", "units", "id_number");
        ncOut_TRACE.addVariableAttribute("fish_identification", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("fish_age", DataType.INT, dim2);
        ncOut_TRACE.addVariableAttribute("fish_age", "units", "days");
        ncOut_TRACE.addVariableAttribute("fish_age", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("Effectif_super_indiv", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("Effectif_super_indiv", "units", "N");
        ncOut_TRACE.addVariableAttribute("Effectif_super_indiv", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("E", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("E", "units", "joules reserve");
        ncOut_TRACE.addVariableAttribute("E", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("V", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("V", "units", "volume structurel");
        ncOut_TRACE.addVariableAttribute("V", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("E_R", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("E_R", "units", "joules reproduction buffer");
        ncOut_TRACE.addVariableAttribute("E_R", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("temperature", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("temperature", "units", "Celsius");
        ncOut_TRACE.addVariableAttribute("temperature", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("f", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("f", "units", "food functional response");
        ncOut_TRACE.addVariableAttribute("f", "_FillValue", Double.NaN);
        
        ncOut_TRACE.addVariable("Nb_eggs_spawned_dt", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("Nb_eggs_spawned_dt", "units", "Number of eggs spawned each day");
        ncOut_TRACE.addVariableAttribute("Nb_eggs_spawned_dt", "_FillValue", Double.NaN);

        ncOut_TRACE.addVariable("E_H", DataType.DOUBLE, dim2);
        ncOut_TRACE.addVariableAttribute("E_H", "units", "joules maturation buffer");
        ncOut_TRACE.addVariableAttribute("E_H", "_FillValue", Double.NaN);


        // create the file
        try {
            ncOut_TRACE.create();
        } catch (IOException e) {
            System.err.println("ERROR creating file " + ncOut_TRACE.getLocation() + "\n" + e);
        }
    }

    static public void writeTRACE_NetCDF(Poissons Po) throws Exception {

        ArrayDouble.D2 longitude_Data = new ArrayDouble.D2(1,drifter_TRACE.getLength());
        ArrayDouble.D2 latitude_Data = new ArrayDouble.D2(1,drifter_TRACE.getLength());
        ArrayDouble.D2 identification_Data = new ArrayDouble.D2(1,drifter_TRACE.getLength());
        ArrayInt.D2 age_Data = new ArrayInt.D2(1, drifter_TRACE.getLength()); // <<<<<<<peut-être viré
        ArrayDouble.D2 Effectif_super_indiv_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength()); // <<<<<<<peut-être viré

// DEB
        ArrayDouble.D2 E_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength());
        ArrayDouble.D2 V_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength());
        ArrayDouble.D2 E_R_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength());
        ArrayDouble.D2 E_H_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength());

        // ENVIRONNEMENT (NON ESSENTIEL)
        ArrayDouble.D2 temperature_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength());
        ArrayDouble.D2 f_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength());


// POUR TESTER LA PONTE :
        ArrayDouble.D2 Nb_eggs_spawned_dt_Data = new ArrayDouble.D2(1, drifter_TRACE.getLength()); // <<<<<<<peut-être viré

       // FILL VALUES
                for (int t = 0; t < drifter_TRACE.getLength(); t++) {
                        longitude_Data.set(0,t, Double.NaN);
                        latitude_Data.set(0,t, Double.NaN);
                        age_Data.set(0,t, -999);
                        Effectif_super_indiv_Data.set(0,t, -999);
                        E_Data.set(0,t, Double.NaN);
                        V_Data.set(0,t, Double.NaN);
                        E_R_Data.set(0,t, Double.NaN);
                        temperature_Data.set(0,t, Double.NaN);
                        f_Data.set(0,t, Double.NaN);
                        identification_Data.set(0,t, -999);
                        Nb_eggs_spawned_dt_Data.set(0,t, Double.NaN);
                }

            int p=0; 
                // On fabrique les data.set :
                longitude_Data.set(0,p, Po.lon);
                latitude_Data.set(0,p, Po.lat);
                identification_Data.set(0,p, Po.id);

                age_Data.set(0, p, Po.age);
                Effectif_super_indiv_Data.set(0, p, Po.S);
                
                E_Data.set(0, p, Po.DEB.E);
                V_Data.set(0, p, Po.DEB.V);
                E_R_Data.set(0, p, Po.DEB.E_R);
                E_H_Data.set(0, p, Po.DEB.E_H);

                temperature_Data.set(0, p, Po.temp);
                f_Data.set(0, p, Po.DEB.f);
                

                Nb_eggs_spawned_dt_Data.set(0, p, Po.Nb_eggs_spawned_dt);

                // On donne l'origine : le pas de temps actuel :
                int[] origin = new int[]{0, 0};
                origin[0] = Poissons.compteur_advec;

                // On ecris les data.set dans le netcdf :
                try {
                    // origin = 0
                    ncOut_TRACE.write("longitude", origin, longitude_Data);
                    ncOut_TRACE.write("latitude", origin, latitude_Data);
                    ncOut_TRACE.write("fish_identification", origin, identification_Data);
                    ncOut_TRACE.write("fish_age", origin, age_Data);
                    ncOut_TRACE.write("Effectif_super_indiv", origin, Effectif_super_indiv_Data);

                    ncOut_TRACE.write("E", origin, E_Data);
                    ncOut_TRACE.write("V", origin, V_Data);
                    ncOut_TRACE.write("E_R", origin, E_R_Data);
                    ncOut_TRACE.write("E_H", origin, E_H_Data);

                    ncOut_TRACE.write("temperature", origin, temperature_Data);
                    ncOut_TRACE.write("f", origin, f_Data);

                    ncOut_TRACE.write("Nb_eggs_spawned_dt", origin, Nb_eggs_spawned_dt_Data);

                } catch (IOException e) {
                    e.printStackTrace();
                }
    //    }// boucle sur pop
    }

    static public void CloseTRACE_NetCDF() throws Exception {

        try {
            ncOut_TRACE.close();
            String strFilePart = ncOut_TRACE.getLocation();
            String strFileBase = strFilePart.substring(0, strFilePart.indexOf(".part"));
            File filePart = new File(strFilePart);
            File fileBase = new File(strFileBase);
            filePart.renameTo(fileBase);
            System.out.println("Closed NetCDF output file.");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
