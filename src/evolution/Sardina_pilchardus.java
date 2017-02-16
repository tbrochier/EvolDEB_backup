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

public class Sardina_pilchardus extends DebLayer{

    public Sardina_pilchardus() {
        // Sardina pilchardus
        // CONDITIONS INITIALES
        E_init = 2.58; // J, Initial Reserve = egg
        V_init = 0.000001; // cm^3, Initial volume --> close to 0, try different initial values
        // taille d'un oeuf = 1 mm3

//Primary parameters        
        // temperature correction
        T_ref = 15 + 273;   //  K, Reference temperature ; 
        T_A = 9800;       //  K, Arrhenius temperature ;
        T_AL = 50000;      //  K, Arrhenius temp for lower boundary
        T_AH = 190000;     //  K, Arrhenius temp for upper boundary

// GAMME MIN-MAX DE TEMP ICI CELLE DES ADULTES INDICATIF,
        //MAIS CHANGE POUR CHAQUE STADE (dans set_stage):
            // Min et Max des corr de flux de temperature pour JUVENILES:
            T_L = 0 + 273;    //  K, lower boundary tolerance range
            T_H = 26 + 273;   //  K, upper boundary tolerance range

       
        // feeding & assimilation
        F_m = 6.51;    // l/d.cm^2, {F_m} max spec searching rate
        kap_X = 0.8;   // -, digestion efficiency of food to reserve
        // -> S. pilachardus : 
        p_Am = 1.677 * 92.51 / 0.3436 * 2.4019;     // J/cm^2/d, maximum surface-specific assimilation rate
        
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
// -> S. pilachardus : 
        E_Hb = 1.372e0;  // J, E_H^b, maturity at birth
        
        E_Hj = 20;       // J, E_H^j, maturity at metamorphosis
        // -> S. pilachardus : 
        E_Hp = 1.928e5;  // J, E_H^p, maturity at puberty

        // param to compute observable quantities
        del_M = 0.1391;    //  -, shape coefficient to convert vol-length to physical length
        d_V = 0.2; 	   // g/cm^3, specific density of structure (dry weight)
        mu_V = 500000;    // J/mol, specific chemical potential of structure
        mu_E = 550000;    // J/mol, specific chemical potential of reserve
        w_V = 23.9;      // g/mol, molecular weight of structure
        w_E = 23.9;     // g/mol, molecular weight of reserve
        c_w = 1 - d_V;   // - , water content (c_w * W_w = total water weight)
        relative_fecondity = 100; // nb oeuf par grammes (A CHERCHER/4 fois inférieur à S. aurita??)

        // compound parameters
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        X_K = Simulation.food_half_saturation; // 0.2 = bon pour le Nano_phyto en surface //(p_Am / (kap_X * F_m))/100;// c'etait 50  // same unit as food, half-saturation coefficient
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        p_Xm = p_Am / kap_X; // J.cm-2.d-1, max surf area specific ingestion rate, here we assume constant assimilation efficiency



        // STAGES-SPECIFIC parameters:
         del_M_egg = del_M; // shape coefficient for EGG
         T_L_egg = T_L;    //  K, lower boundary tolerance range
         T_H_egg = T_L;  //  K, upper boundary tolerance range

            del_M_YOLK_SAC_LARVA = del_M; // shape coefficient for YOLK_SAC_LARVA
            T_L_YOLK_SAC_LARVA = T_L;//0 + 273;    //  K, lower boundary tolerance range
            T_H_YOLK_SAC_LARVA = T_H;//27 + 273;   //  K, upper boundary tolerance range
                
            del_M_FEEDING_LARVA = del_M;// shape coefficient for FEEDING_LARVA
            T_L_FEEDING_LARVA = T_L;    //  K, lower boundary tolerance range
            T_H_FEEDING_LARVA = T_H;//27 + 273;   //  K, upper boundary tolerance range
               
            del_M_JUVENILE = del_M;// shape coefficient for JUVENILE
            T_L_JUVENILE = T_L;    //  K, lower boundary tolerance range
            T_H_JUVENILE = T_H;//27 + 273;   //  K, upper boundary tolerance range
            
            del_M_ADULT = del_M;// shape coefficient for ADULT
            T_L_ADULT = T_L;   //  K, lower boundary tolerance range
            T_H_ADULT = T_H;//27 + 273;   //  K, upper boundary tolerance range
 
    
    
    
    }
        
}

