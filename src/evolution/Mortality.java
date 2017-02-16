/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package evolution;

/**
 *
 * @author timbrochier
 */
public class Mortality {

    static float Mmax = 0.5f; // Mortalite naturelle max par jour
    static float Mmin = 0.000f; // Mortalite naturelle min par jour (pr Lmax 44,8 = 0.001f

    static float Mmax_compet = 0.1f; //Mortalite de competition max par jour
//static float param_ajsutement_predation = 0.7f; // plus il est grand plus la mortalite diminue
    static double M_predation, M_competition, M_senescence, M_peche;
    static int Biom_landings;

    public static double calcul_mortality(float fish_length, double fish_weigth, int fish_stade, int fish_age, double effectif, float biomasse_autres_SI_1, float bathy, float dt_jour, int zone) {

        // Taille des oeufs :
        float L_init = 0.1f; //cm

// 1) COMPETITION
        //M_competition = mort_competition_function(effectif, fish_weigth, Carcapa, biomasse_autres_SI_1);    
//        if (bathy>-200){            
//        M_competition = mort_competition_function(zone, Carcapa);
//        }else{
//           M_competition = 0;
           // (considerons que au large il n'y a pas de pb de competition por l'espace)
//        }

//if   (M_competition > 0){
        //         System.out.println(" COMPETITION! M_competition = " + M_competition + " ; BIOMASSE de ce SI = " + effectif*fish_weigth/1e6 + " tonnes" );    
        //      System.out.println(" COMPETITION! M_competition = " + M_competition + " ; biomasse_autres_SI_1 = " + biomasse_autres_SI_1 + " ; Carcapa = " +  Carcapa +  " effectif de ce SI = "  + effectif  +" ; fish_weigth = " +  fish_weigth );    
//} 
// densite_biomasse_locale= somme des biomasses de ce bancs et des autres dans la maille de grille, divise par la surface de la maille en km2
// Carcapa = somme du carbonne disponible par km2 sur les 20 premiers metres de la colonne d'eau
        double biomasse_cell_TOTAL = biomasse_autres_SI_1 + fish_weigth * effectif;
        /*
         if  (biomasse_cell_TOTAL > 1.2*Carcapa){
         // (// 1.2 : on autorise un dépassement de 20% de la carcapa, sinon on risque
         // de faire des calculs de competition sur des effectif de surfplus très faible,
         // avec meme le risque de division par zero si round(effectif_surplus)=0

         System.out.println(" COMPETITION! biomasse TOTALE DANS LA CELLULE = " + biomasse_cell_TOTAL/1000000 + " tonnes ; Carcapa = " + Carcapa/1000000 +" tonnes ");
         //    double capamax; // effectif max pour une cohorte d'indiv de ce poids (fish_weigth)
         double surplus_biomasse=biomasse_cell_TOTAL-Carcapa; // biomasse en surplus (qui depasse la Carcapa)
    
         //double effectif_pour_cacapa = Carcapa/fish_weigth;
         double surplus_effectif=Math.round(surplus_biomasse/fish_weigth); // effectif en surplus (qui depasse la capamax)

         //    M_competition = mort_competition_function(surplus_effectif);
         //    effectif = effectif_pour_cacapa;//Math.max(effectif_pour_cacapa,effectif-surplus_effectif);// + ( Math.round(surplus_effectif*(1-Mmax_compet*M_competition*dt_jour)));
    
         //    effectif = Math.max(0,effectif-surplus_effectif) + ((int) Math.round(surplus_effectif*(1-Mmax_compet*M_competition*dt_jour)));
         }
         */

// 2) MORTALITE PAR PREDATION (SELON LA TAILLE):
// Correctif bathy : pour représenter le fait qu'il y a plus de prédateurs en haute mer, et + gros,
// on considere un facteur correctif à la taille du poisson avant calcul de la mortalite par predation :
        float facteur_correctif_mort_bathy = Math.max(0.1f, 1 - (float) Kinesis.predation_bathy(bathy));
        float fish_length_corr = facteur_correctif_mort_bathy * fish_length;
        M_predation = mort_size_function(fish_length_corr, fish_stade, fish_age);

        if (bathy < -20000) {
            System.out.println("bathy = " + bathy + " --> facteur_correctif_mort_bathy = " + facteur_correctif_mort_bathy + " ; M_predation = " + M_predation);
            System.out.println("fish_length_corr = " + fish_length_corr + " ; fish_stade = " + fish_stade + " ; dt_jour = " + dt_jour);
            System.out.println("..   ");
        }

        double beta = 2e-13; // 2e-13 jusqu'au 27 juin (20e-13 ensuite, jusqu'au 5 octobre 2015) 
//double M_sen = Math.min(1, beta*fish_age*fish_age*fish_age);
//M_senescence = Math.max(M_sen, Mmin);
        M_senescence = Math.min(1, beta * fish_age * fish_age * fish_age);
//if (bathy<-2000){
//System.out.println("bathy = " + bathy + " --> facteur_correctif_mort_bathy = " + facteur_correctif_mort_bathy);
//}

// 3) Mortalite de peche
// sur tout individu > 10cm qui se trouve sur le plateau continental
        M_peche = mort_peche_function(fish_length); //, bathy_actuelle)*dt_jour;

// competition : 
        effectif = Math.floor(effectif * (1 - M_competition * dt_jour));
// Predation : 
        effectif = Math.floor(effectif * (1 - M_predation * dt_jour));
// Senescence : 
        effectif = Math.floor(effectif * (1 - M_senescence * dt_jour));
// Peche : 
        Biom_landings = (int) Math.round(effectif * M_peche * fish_weigth); // en grammes
        effectif = Math.floor(effectif * (1 - M_peche * dt_jour));

//effectif =  Math.round( effectif*(1 - Mmax*(M_predation + M_senescence) - M_peche));
//            System.out.println(" nouvel_effectif = " + nouvel_effectif);
// System.out.println(" bathy = " + bathy + " ; facteur_correctif_mort_bathy = " + facteur_correctif_mort_bathy);
        // On prend un L_init qui correspond a la metamorphose = 60 jours
        return effectif;
    }

    static double mort_peche_function(float fish_length) {
    // il faudrait faire dependant du stade (capturabilite)
        // + il faut enregistrer les tonnes de captunes par zone
        double M_pech;

        if (fish_length > Simulation.Min_fishing_size) {
            M_pech = Simulation.F_annuel / (float) Simulation.nb_jours_par_ans;
        } else {
            M_pech = 0;
        }
        return M_pech;
    }

    ;

static double mort_size_function(float fish_length, int stage, int age) {
        double M_pred;

    //float epsilon = L_init/10;
//    double M = (L_init-epsilon)/fish_length; // mortalite indicative, a multiplier par Mmax pou l'avoir par jour
        //double M = Math.exp((L_init-epsilon)/fish_length) / Math.exp((L_init-epsilon)/L_init); // mortalite indicative, a multiplier par Mmax pou l'avoir par jour
    // Okunishi et al 2012 // voir aussi refs dans Yi Xu et al 2013
        if (stage == 0) { // Eggs
            M_pred = 0.57;//0.9 (pour Lmax 44,8); //0.57;
        } else if (stage == 1) { // Yolk sac larvae
            M_pred = 0.3; //0.3;
        } else { // a partir du stade de feeding larvae, mortalite dependante de la taille : 
            M_pred = 0.189 * Math.exp(-fish_length / 2.468);
//M_pred = 0.9*Math.exp(-fish_length/3); //25 juin 2014
        }

        M_pred = Math.max(Mmin, M_pred);
        return M_pred;
    }
    /*
     static double mort_competition_function(double effectif, double fish_weigth, double Carcapa, double biomasse_autres_SI){
     double MM;

     double carcapa_locale_restante = Carcapa - biomasse_autres_SI;    


     if (carcapa_locale_restante>0){
     double effectif_pour_cacapa_restante = carcapa_locale_restante/fish_weigth;

     if (effectif>effectif_pour_cacapa_restante){ // COMPETITION

     // System.out.println(" COMPETITION! effectif = " + effectif + "  ; effectif_pour_cacapa_restante = " + effectif_pour_cacapa_restante + "  ; Carcapa = " + Carcapa/1000000 + " tonnes ; fish_weigth = " + fish_weigth + " grammes");
     MM = 1 - effectif_pour_cacapa_restante/effectif;

     }else{ // PAS DE COMPETITION
     MM=0;
     }
     }
     else{ // si deja plus de place ici : 
     MM = Mmax_compet;
     }
     MM = Math.min(Mmax_compet, MM);
     // MORTALITE QUADRATIQUE :
     //double alpha = 1/(capamax+1);
     //double MM = Math.min(1,alpha*effectif);

     // MORTALITE EN SERIE f(n) = 1 - 1/(2^n)
     //  int n = (int) Math.round(surplus_effectif/1000);
     //     MM = 1 - 1/(Math.pow(2, n)) ;

     return MM;
     }

     static double mort_competition_function(double surplus_effectif){
     double MM;

     // MORTALITE EN SERIE f(n) = 1 - 1/(2^n)
     int n = (int) Math.round(surplus_effectif/1000);

     MM = 1 - 1/(Math.pow(2, n)) ;

     return MM;
     }

     */

    static double mort_competition_function(int zone, double carcapa_zone_locale) {
        double MM;
        if (Population.Biom_zone[zone] > carcapa_zone_locale) {
            MM = Population.Biom_zone[zone] / carcapa_zone_locale - 1;
            System.out.println(" COMPETITION! Population.Biom_zone[zone] = " + (Population.Biom_zone[zone]/1e6) + " tonnes ; carcapa_zone_locale = " + (carcapa_zone_locale/1e6) + " tonnes");
        } else {
            MM = 0;
        }
        return MM;
    }
}
