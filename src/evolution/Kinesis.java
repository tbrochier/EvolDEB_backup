package evolution;

import java.util.Random;
import static java.lang.Double.*;
import evolution.DebLayer;
/**
 *
 * @author timbrochier
 * Adapted from Kinesis fish movement described in Watkins & Rose 2013
 */
public class Kinesis {

    static double[] V_move = new double[3]; // V_move[0] = Vx ; V_move[1] = Vy; // Deplacement en x et y
    static double M; // mortalité : doit etre static pour pouvoir etre appele dans poisson
    static double M_Dens; // index de mortalite de competition : public car enregistrer dans sorties
    static double Mmax = 1;//0.5 / (24*60 / Simulation.dt_advec); // 0.5 par jour
//        %1.2/(24/Simulation.dt_advec); // C'est la valeur max de la mortalité qu'on va appliquer aux super-individus a chaque pas de temps dt_advec
// Houde, 2008 : Mortalite max par jour pour les larves est de 1,2 (S = S*exp(-1.2) ==> 70% de mort)
    static double Gmax = (0.3 / 24*60) * Simulation.dt_advec;  //grammes/heure * dt_advec (en heure)   Croissance max en gramme par pas de temps...
// AU PIF : mettons 0.3 g par jours !!!
    static double mDDmax_factor = 10;//
    static double Lmax = 25; //cm - La taille max des poissons, pour definir la mortalite size-dependant
    static double L_init = 0.3; // La taille min des poisson (3mm? --> CHECK)
// POUR LES 5 PARAMETRES SUIVANT : WATKINS ET ROSE 2013 - VALEURS ALGO GENETIQUE!!
    static double delt = Simulation.delta_kinesis; //0.2;//0.5;//0; //pondere growth and local_mortality dans Q : 0-> tout sur Growth; 1-> tout sur mortalite
    static double H1 = 0.6; // hauteur de la fonction d'inertie
    static double H2 = 1; // hauteur de la fonction de random
    static double H3 = 0.5; // hauteur de la fonction de extended Kinesis
    static double Qopt = 1; // valeur de optimal habitat
    static double sigmaQ = 0.1299; //0.2; //Plage de changement du comportement (inertie/rendom). sigmaQ petit --> changement sharp
    static Random myRandom = new Random();
    static int[] signe = new int[]{-1, 1};

    static double grad_temp_ok = 0.1; // degre °C par km. Des gradients plus fort auront un effet de repousser les sardinelles

    public static double Habitat_quality(Poissons Po){
            /*
            double effectif, float Densite_biomass, double CarCapa, float fish_size,
            double temperature, double grad_temp, double f, double  T_opt, double sigma_t, double lon, double lat, int bathy) {
*/
        double Q;// Qualite de l'habitat

        //double G = delta_weigth; // PAS BON CAR weigth diminue quand ponte,
        // même si l'environnement est bon..

        // 16 fev 2016 : dans indice_temperature_bouffe je passe Po.SST au lieu de Po.temp
        double I_tf = indice_temperature_bouffe(Po.SST, Po.grad_temp, Po.f_reponse_fonctionnel,  Po.T_opt, Po.sigma_t);
        int zone = (int) Math.round(Po.lat - Dataset_EVOL.lat_min);

        M = local_mortality(Po.S, Po.weight, Po.biomasse_locale_autres_SI, Po.Body_length, Po.bathyActuelle, zone);
    
        double G_prime = I_tf;//growth_cue(G, Gmax);
        double M_prime = mortality_cue(M);
//System.out.println("Po.bathyActuelle = " + Po.bathyActuelle + " M_prime = " + M_prime);

        // !!!!!! POUR TEST TESTESTESTEST
        // TEST1 = Q = f;//delt + (1-delt)*G_prime - delt*M_prime;
        // il y avait pb dans ces test : on prenait la bouffe a z=-1 (=0 = fond)
        // et en plus f= 0 pour les 1er pas de temsp (larve)
        // TEST2 :bouffe a z=31 (surface)
        // soit G_prime = f (pour la phase de validation TEST_KINESIS)
        G_prime = I_tf;

        Q = delt + (1-delt)*G_prime- delt*M_prime;

        /*
if (Densite_biomass/1000000>5){
System.out.println("DEPLACEMENT INFLUENCE PAR DENSITE EN CET ENDROIT LA; Dens_biom = " + Densite_biomass/1000000 + " tonnes par km2 " );
System.out.println(" mDD = " + M_prime + "; f =  "+ f +" ; Q = " + Q);
        }
*/
//System.out.println("Po.bathyActuelle = " + Po.bathyActuelle + " Q = " + Q);
        
        // DEBUG
        if (isNaN(Q)) {
            System.out.println("BUG KINESIS Q - Densite_biomass =  " + Po.biomasse_locale_autres_SI
                    + " , CarCapa = " + Po.CarCapa + " , fish_size = " + Po.Body_length
                    + " , G_prime = " + G_prime + " , M_prime = " + M_prime);
            System.out.println("temperature = " + Po.temp + " ; grad_temp = " + Po.grad_temp +
                    " ; f = " + Po.f_reponse_fonctionnel + " ; T_opt = " + Po.T_opt + "sigma_t = " + Po.sigma_t);
        }
        return Q; // Qualite de l'habitat
    }

    static double growth_cue(double G, double Gmax) {
        double G_pr;
//    G = CHANTIER EN COURS
        G_pr = G / Gmax;
        return G_pr;
    }

    static double mortality_cue(double M) {
        double M_pr;
        M_pr = (1 - Math.exp(-M)) / (1 - Math.exp(-Mmax));
//        System.out.println("M_locale = " + M + " ; mortality_cue = " + M_pr);                
        return M_pr;
    }

    /*
    static double growth_change(double densE, double densE_old){
    double G_t = densE - densE_old;
    return G_t;
    }
     *
     */
    
    static double local_mortality(double effectif, double fish_weigth, double Dens_biom_autres, float fish_size, float bathy, int zone) {

        // 0 = pas de mortalite
        // 1 = mortalite maximale (100% de morts!)
        
        // TOUT CE QUI PEUT IMPACTER LA MORTALITE LOCALEMENT,
        // (--> autre que par le bilan energetique individuel <--)

        // 1) La Carrying capacity : density dependance
//if (Population.Biom_zone[zone]>CarCapa){
//      System.out.println("Zone très dense : zone " + zone );                
// System.out.println("SI se croisent... biom des autres = " + Dens_biom_autres/1E6 +" tonnes");                
//}

//// EN CHANTIER : 
//        if (Simulation.FLAG_KinesisDD && Dens_biom_autres>0){
////            M_Dens = mort_DD(effectif, fish_weigth, Dens_biom_autres, CarCapa);            
//            //disons
//            M_Dens=1;
//        }else{M_Dens=0;}
//        if (M_Dens>0.5){
//     //  System.out.println("On fuis cet endroit car il y a déja trop de monde!!");            
//        }
//// M_DD = 0 si Dens_biom<CarCapa
        // a partir de Dens_biom=CarCapa, augmente jusqu'a
        // M_DD = 1 si Dens_biom=>2*CarCapa

        // 2) LA PECHE
        // M_peche

        // 3) PREDATION
      //  double M_peche = mort_Pred_peche(bathy);
        // A REVOIR : PLUTOT DISTANCE A LA COTE QUE BATHY (QTTE D'ESSENCE POUR
        // PECHE ARTISANALE)
        // 4) L'instinct de Patrice ("Indice bathymetrique")

double  M_ibat = predation_bathy(bathy);
   //   System.out.println("Bathy = " + bathy + " ; M_ibat =  " + M_ibat );
//      double M_local = M_Dens*M_peche*M_Ibat;//*M_peche*M_pred M_DD*
        
// Le module de densite ne fonctionne pas: (--> seule mortalite locale = M_ibat)
M_Dens = 0;
double M_local = M_ibat + M_Dens - M_ibat*M_Dens;//*M_peche*M_pred M_DD*
        double M_tot = Math.min(Mmax , M_local);
    //  System.out.println("M_tot =  " + M_tot );
        return M_tot;
    }

    static double mort_DD(double effectif, double fish_weigth, double Dens_biom, double CarCapa) {
        double M_norm;

        // PRENDRE IC LA MEME MORTALITE DE COMPETITION QUE DANS Mortality.mort_competition_function( effectif,  fish_weigth,  Carcapa, biomasse_locale_autres_SI)
        M_norm = Mortality.M_competition;
                //Mortality.mort_competition_function(effectif, fish_weigth, CarCapa, Dens_biom);
//        ATTENTION IL FAUT AVOIR APPELER mort_competition_function avant la Kinesis
        // ancien : 
        /*
        double r = Dens_biom/CarCapa;
        double sigma_r = 0.3;
        mDD = Math.exp(0.5*Math.pow(r/sigma_r,2));
        M_norm = Math.min((mDD / mDDmax_factor), 1);
//        int CarCapa_cor = (int) (CarCapa/mDDmax_factor);
        // CarCapa_cor = c'est le "crowd" percut par les sardinelles
  //      if (Dens_biom < CarCapa_cor) {
    //        mDD = 0;
      //  } else {
        //    mDD = (Dens_biom / CarCapa_cor) - 1;
  //
        // M_norm = entre 0 et 1
*/
        return M_norm;
    }

    static double mort_Pred_peche(float bathy) {
        double M_pp=0;
        if (bathy>1000){
           M_pp = 0.5;
        }
        return M_pp;
    }

    public static double[] move(double[] previous_V, float fish_length, double Q, double Qold) {
        double deltaQ;
 // ATTENTION, il faut bien que "previous_V" soit en [m/s]
        double phi = Simulation.Vitesse_max_Bls * fish_length / 100;// C'est la vitesse de nage max en metre par seconde

        if (Q >= Qopt) {
            deltaQ = 0;
        } else {
            deltaQ = Q - Qopt;
        }

// "Habitat condition" : varie entre 0 (habitat nul) et 1 (habitat optimal)
        double Ih = Math.exp(-0.5 * Math.pow(deltaQ / sigmaQ, 2));

        if (Qold >= Q) {
            // KINESIS NORMALE (WATKINS ET ROSE 2012)
            double fx = inertie(previous_V[0], Ih);
            double fy = inertie(previous_V[1], Ih);
            double fz = inertie(previous_V[2], Ih);

            double gx = random(phi, Ih); // [m/s]
            double gy = random(phi, Ih); // [m/s]
            double gz = random(phi, Ih); // [m/s]

// Prenons la moyennes des terme d'inertie et de random$
            // car [bornes inertie] = [-phi phi]
            // et  [bornes random]  = [-phi phi]
            // donc il gaut prendre la moyenne et pas la somme pour avoir une 
            // vitesse bornée sur [-phi phi] sur chaque direction
            // (approximation de la distance reelle = sqrt(Vx^2 + Vy^2) 
            
            double Vx = (fx + gx)/2;  // [m/s]
            double Vy = (fy + gy)/2; // [m/s]
            double Vz = (fz + gz)/2; // [m/s]

            V_move[0] = Vx;
            V_move[1] = Vy;
            V_move[2] = Vz;

            // DEBUG
        if ( (Math.abs(V_move[0])> phi) || (Math.abs(V_move[1])> phi)){                            
            System.out.println("fx = " + fx + " ; fy = " + fy + " ; phi = " + phi);
            System.out.println("gx = " + gx + " ; gy = " + gy);
        }
        
            if (isNaN(Vx)) {
                System.out.println("BUG KINESIS 1 - fish_length =  " + fish_length + " , Q = " + Q + " , Qold = " + Qold);
            }

        } else if (Qold < Q) {
// EXTENTED KINESIS :   OKUNISHI et al. 2012
            // SI L'HABITAT S'AMELIORE ON GARDE LA MEME DIRECTION
            V_move[0] = previous_V[0] - (1 - (Math.abs(previous_V[0]) / phi) * H3) * previous_V[0] * Ih;
            V_move[1] = previous_V[1] - (1 - (Math.abs(previous_V[1]) / phi) * H3) * previous_V[1] * Ih;
            //V_move[2] = previous_V[2] - (1 - (Math.abs(previous_V[2]) / phi) * H3) * previous_V[2] * Ih;

        // CORRECTION 23 juin : 
    /*        
            if (previous_V[0]>0){
            V_move[0] = previous_V[0] - (1 - phi/(Math.abs(previous_V[0])) * H3) * previous_V[0] * Ih;                
            }else{V_move[0]=0;}
            if (previous_V[1]>0){
            V_move[1] = previous_V[1] - (1 - phi/(Math.abs(previous_V[1])) * H3) * previous_V[1] * Ih;
            }else{V_move[1]=0;}
*/        }
            

        // DEBUG
/* PAS BON (phi tjr >0)
        if (Math.abs(V_move[0])> phi){
            V_move[0] = phi;}
        if (Math.abs(V_move[1])> phi){
            V_move[1] = phi;}
  */          
        if ( (Math.abs(V_move[0])> 2*phi) || (Math.abs(V_move[1])> 2*phi)){
        V_move[0] = V_move[0]/10;
        V_move[1] = V_move[1]/10;
        
//    System.out.println("BUG KINESIS 2 - fish_length =  " + fish_length + " , Q = " + Q + " , Qold = " + Qold);
  //          System.out.println("V_move[0] = " + V_move[0] + " , V_move[1] = " + V_move[1] + " ; phi = " + phi);
    //        System.out.println("previous_V[0] = " + previous_V[0] + " ; previous_V[1] = " + previous_V[1]);
      //      System.out.println("Ih = " + Ih);
        }
        return V_move;
    }

    static double inertie(double previous_V, double Ih) {
        double f = previous_V * H1 * Ih;
        return f;
    }

static double indice_temperature_bouffe(double temperature, double grad_temp, double f, double  T_opt, double sigma_t){

    // Flux de temperature AREHENIUS modifie :
    // modif : si flux > 1, on considere negatif car metabolisme
    // doit tourner trop vite
//    double It = DebLayer.getcT(temperature);
//    if (It>1){ It = 1 - (It-1);}
    // PRENONS PLUTÖT UNE GAUSSIENNE CENTREE AUTOUR DE LA TEMPERATURE DE NAISSANCE :
//double  T_opt; // Temperature cible par ce poisson
//double sigma_t; // largeur de la gaussienne autour de T_opt

double It = Math.exp( - Math.pow(temperature-T_opt,2) / 2*Math.pow(sigma_t,2) );

if (Simulation.Flag_pref_thermique_nage==false){
    It = 1;
}
if (Simulation.Flag_pref_bouffe_nage==false){
    f = 1;
}
    // aussi, les sardinelles detestent les changements brusques de temperature
    // (voir ref Sardinelles mortes en Grece)
    // si grad_temp est + grand que grad_temp_ok, alors ça diminuera l'indice
    // d'habitat et ça fera tourner la sardinelle
    //    double I_grad_temp = Math.min(grad_temp_ok/Math.abs(grad_temp), 1);

    double I = f*It;//*I_grad_temp;
    return I;
}


    static double random(double phi, double Ih) {

        // --> 1 body_length per seconde => convertir fish_length en metre

        double moy = Math.sqrt(phi * phi / 2);
        double std = phi / 2;

        int rdm = myRandom.nextInt(signe.length);

        double epsilon = getGaussian(moy, std) * signe[rdm];
        double g = epsilon * (1 - H2 * Ih);

        // 18 mai 2015 : 
        // ajoutons une borne de securite : en effet, meme si c'est tres rare, 
        // epsilon peut prendre des valeur trop grande
        // assurons nous donc que g < phi (phi = vitesse max)
        if (Math.abs(g)> phi){ g = phi * signe[rdm];}
        
        /*            boolean Iwant = true;
        while (Iwant){
        System.out.println(epsilon);
        rdm = myRandom.nextInt(signe.length);
        epsilon = getGaussian(moy, std) * signe[rdm];            }
         */
        return g;
    }

    /**
    Generate pseudo-random floating point values, with an
    approximately Gaussian (normal) distribution.

    Many physical measurements have an approximately Gaussian
    distribution; this provides a way of simulating such values.
     */
    public static double getGaussian(double aMean, double aVariance) {
        Random fRandom = new Random();
        return aMean + fRandom.nextGaussian() * aVariance;
    }



public static double predation_bathy(float bathy){

    int plateau = -200;
    int bat0 = -2000;

    double M_ibat_lineaire, M_ibat;
    if (bathy>plateau){
        M_ibat_lineaire = 0;
    }else {
    M_ibat_lineaire =  Math.min(1, ((float) bathy)/((float) bat0) - ((float) plateau)/((float) bat0));
    }
double sigma_M_Ib_lineaire = 0.3;
M_ibat = Math.exp(  -0.5* Math.pow((1-M_ibat_lineaire)/sigma_M_Ib_lineaire , 2) );


//System.out.println("bathy =  " + bathy + " M_ibat =  " + M_ibat);

return M_ibat;
}

}



/*
    static double mort_size(double size_fish) {

        double l = Math.max(L_init, size_fish);
//    double ML = (Lmax - l) / (Lmax - L_init); (Watkins et Rose)
        double ML = (1 / l) / (1 / L_init); //(moi meme : mortalite diminue de façon hyperbolique avec la taille)
//    System.out.println("ML =  " + ML + " size_fish =  " + size_fish);

        return ML;
    }
*/