/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package evolution;

/**
 *
 * @author timbrochier
 */

public class Sardinella_aurita extends DebLayer{

    public Sardinella_aurita() {
        // CONDITIONS INITIALES
        E_init = 2.58; // J, Initial Reserve = egg
        V_init = 0.000001; // cm^3, Initial volume --> close to 0, try different initial values
        // taille d'un oeuf = 1 mm3

        // STAGES-SPECIFIC parameters:
         del_M_egg = 0.1391; // shape coefficient for EGG
         T_L_egg = 15 + 273;    //  K, lower boundary tolerance range
         T_H_egg = 27 + 273;  //  K, upper boundary tolerance range

            del_M_YOLK_SAC_LARVA = 0.1391; // shape coefficient for YOLK_SAC_LARVA
            T_L_YOLK_SAC_LARVA = 15 + 273;//0 + 273;    //  K, lower boundary tolerance range
            T_H_YOLK_SAC_LARVA = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range
                
            del_M_FEEDING_LARVA = 0.1391;// shape coefficient for FEEDING_LARVA
            T_L_FEEDING_LARVA = 15 + 273;    //  K, lower boundary tolerance range
            T_H_FEEDING_LARVA = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range
               
            del_M_JUVENILE = 0.1391;// shape coefficient for JUVENILE
            T_L_JUVENILE = 15 + 273;    //  K, lower boundary tolerance range
            T_H_JUVENILE = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range
            
            del_M_ADULT = 0.1391;// shape coefficient for ADULT
            T_L_ADULT = 10 + 273;   //  K, lower boundary tolerance range
            T_H_ADULT = 29 + 273;//27 + 273;   //  K, upper boundary tolerance range
 
//Primary parameters        
        // temperature correction
        T_ref = 20 + 273;   //  K, Reference temperature ; 
        T_A = 5000;//9800;       //  K, Arrhenius temperature ;
        T_AL = 50000;      //  K, Arrhenius temp for lower boundary
        T_AH = 190000;     //  K, Arrhenius temp for upper boundary

// GAMME MIN-MAX DE TEMP ICI CELLE DES ADULTES INDICATIF,
        //MAIS CHANGE POUR CHAQUE STADE (dans set_stage):
            // Min et Max des corr de flux de temperature pour JUVENILES:
            T_L = 15 + 273;//0 + 273;    //  K, lower boundary tolerance range
            T_H = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range

        
        // Relation d'echelle entre S. pilchardus et S. aurita
        // (en considerant un seul et meme shape coef pour les deux : 
        // del_M = 0.1391;)
        // --> on a donc le zoom factor suivant z =  Lm(sardinelle) / Lm (sardine) = 1.5471
        

        // feeding & assimilation
        F_m = 6.51;    // l/d.cm^2, {F_m} max spec searching rate
        kap_X = 0.8;   // -, digestion efficiency of food to reserve
        // -> S. pilachardus : p_Am = 1.677 * 92.51 / 0.3436 * 2.4019;     // J/cm^2/d, maximum surface-specific assimilation rate
        // -> S. aurita :
        // pour Longueur max de 44,8cm : 
//        p_Am = 1.6778e+03;     // J/cm^2/d, maximum surface-specific assimilation rate        
        // pour Longueur max de 40cm :         
        p_Am = 1.4980e+03; // J/cm^2/d, maximum surface-specific assimilation rate

// mobilisation, maintenance, growth & reproduction
        v = 0.1379 * 2.4019;  // cm/d, energy conductance
        Kappa = 0.3436;       // -, allocation fraction to soma = growth + somatic maintenance
        kap_R = 0.95;         // -, reproduction efficiency
        p_M = 92.51;          // J/d.cm^3, [p_M], vol-specific somatic maintenance
        p_T = 0;              // J/d.cm^2, {p_T}, surface-specific som maintenance
        k_J = 0.002;          // 1/d, maturity maint rate coefficient
        E_G = 4767;          // J/cm^3, [E_G], spec cost for structure

        // life stages: E_H is the cumulated energy from reserve invested in maturation
        E_Hh = 1;        // J, E_H^h, maturity at hatching
// -> S. pilachardus : E_Hb = 1.372e0;  // J, E_H^b, maturity at birth
        //-> S. aurita :
        // pour Longueur max de 44,8cm : 
       // E_Hb = 1.2; //5.0805;  // J, E_H^b, maturity at birth / S. aurita
        // pour Longueur max de 40cm : 
        E_Hb = 1.2; //3.6159;  // J, E_H^b, maturity at birth / S. aurita
        
        
        E_Hj = 20;       // J, E_H^j, maturity at metamorphosis
        // -> S. pilachardus : E_Hp = 1.928e5;  // J, E_H^p, maturity at puberty
        // pour Longueur max de 44,8cm : 
//        E_Hp = 7.1394e5;//1.6000e+05;//7.1394e5;  // J, E_H^p, maturity at puberty
        // Changement : de 7.1394e5 a 1.6000e+05 le 19 juin 2015
        // 21 juin : retour a 7.1394e5 car sinon on a des taille max <20cm!
        // pour Longueur max de 40cm : 
        E_Hp = 5.0813e+05;  // J, E_H^p, maturity at puberty

        // essais E_Hp = 1.928e6 simu 39 --> pas de maturation avant 2 ans (deux ans de run)
        // essais E_Hp = 1.928e5 simu 39 --> pas de maturation avant 6 mois

        // param to compute observable quantities
        del_M = 0.1391;    //  -, shape coefficient to convert vol-length to physical length
        d_V = 0.2; 	   // g/cm^3, specific density of structure (dry weight)
        mu_V = 500000;    // J/mol, specific chemical potential of structure
        mu_E = 550000;    // J/mol, specific chemical potential of reserve
        w_V = 23.9;      // g/mol, molecular weight of structure
        w_E = 23.9;     // g/mol, molecular weight of reserve
        c_w = 1 - d_V;   // - , water content (c_w * W_w = total water weight)
        relative_fecondity = 400; // nb oeuf par grammes de femelle These Freon

        // compound parameters
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        X_K = Simulation.food_half_saturation; // 0.2 = bon pour le Nano_phyto en surface //(p_Am / (kap_X * F_m))/100;// c'etait 50  // same unit as food, half-saturation coefficient
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        p_Xm = p_Am / kap_X; // J.cm-2.d-1, max surf area specific ingestion rate, here we assume constant assimilation efficiency

    }
        
}

