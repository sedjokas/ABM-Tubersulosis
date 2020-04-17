/***
 *  Name: TBmicro
 *  Author: Selain Kasereka
 *  E-mail: selain.kasereka@unikin.ac.cd
 *  Supervision: Emile-Franc Doungmo Goufo and Ho Tuong Vinh
 *  Description: Agent-Based Simulation of TB spreding, case study of DRC. This model is based on a compartmental model proposed by Selain Kasereka
 *  in an article concerning mathematical modeling and simulation of TB in RD Congo.
***/

model TBmicro1

global {
	//Number of susceptible individuals
	int S_people <- 700; //Number of susceptible individuals
	int I_people <- 20;  //Number of infectious individuals
	int Le_people <- 0; //Number of Latent Early individuals
	int Lf_people <- 0; //Number of Latent Late individuals
	int R1_people <- 0; //Number of Recovered Spontaneously individuals
	int R2_people <- 0; //Number of Recovered after traitment
	int T_people <- 0; //Number of Transfered individuals
	int K_people <- 0; //Number of Lost of Sight individuals
	
	//Radius of contamination
	int range;
	
	//Step of simulation
	//float step<- 1#year;
	
	//Total number of individuals
	int N <- S_people + I_people + Le_people + Lf_people + R1_people + R2_people + T_people + K_people;
	
	
	//My parameter
	float Gamma ; // Rate of recruitment in compartment S 
	float mu1 ; //rate of mortality not related to TB infection 
	float mu2 ; //rate of mortality related to TB infection  
	float beta ; //rate of transferred people in a other hospital
	float gamma ; //rate of recovered after treatment process
	float sigma ; //rate of spontaneously recovered
	float lambda ; // rate of transmission
	float q ; // rate of progression to active TB (Le ->I)
	float r ; //rate of re-infection from R1 and R2 to Le
	float r1 ; //rate of re-infection from I to Le
	float r2 ; //rate of re-infection from T to Lf
	float r3 ; //rate of re-infection from K to Le
	float g1 ; //rate of recovered after treatment process from Le to R1
	float g2 ; // rate of spontaneously recovered from Le to R2
	float k1 ; // rate of recovered after treatment process from Lf to R1 
	float k2 ; //rate of spontaneously recovered from Lf to R2
	float h ; //rate of progression of TB infection, from Le to Lf
	float w ; // rate of slow progression of TB infection, from Lf to I 
	float v ; // rate of infected who interrupt their treatment 
	

	init {
		//Create the number of hosts susceptibles
		create Host number: floor(S_people) {
			is_susceptible <- true;
			is_infected <- false;
			is_latent_e <- false;
			is_latent_f <- false;
			is_recovered_1 <- false;
			is_recovered_2<- false;
			is_transfered <- false;
			is_lostsight <- false;
			color <- #green;
			//is_infected_sinse <- 0.0;
		}
		
//Create the numberof Infected individuals
		create Host number: I_people {
			is_susceptible <- false;
			is_infected <- true;
			is_latent_e <- false;
			is_latent_f <- false;
			is_recovered_1 <- false;
			is_recovered_2<- false;
			is_transfered <- false;
			is_lostsight <- false;
			color <- #red;
		}
		
		//Create the number of Latent Early hosts 
		create Host number: Le_people {
			is_susceptible <- false;
			is_infected <- false;
			is_latent_e <- true;
			is_latent_f <- false;
			is_recovered_1 <- false;
			is_recovered_2<- false;
			is_transfered <- false;
			is_lostsight <- false;
			color <- #orange;
			
		}
		
		//Create the number of Latent Late hosts 
		create Host number: Lf_people {
			is_susceptible <- false;
			is_infected <- false;
			is_latent_e <- false;
			is_latent_f <- true;
			is_recovered_1 <- false;
			is_recovered_2<- false;
			is_transfered <- false;
			is_lostsight <- false;
			color <- #yellow;
		}
		
		//Create the number of Recovered after traitment hosts 
		create Host number: R1_people {
			is_susceptible <- false;
			is_infected <- false;
			is_latent_e <- false;
			is_latent_f <- false;
			is_recovered_1 <- true;
			is_recovered_2<- false;
			is_transfered <- false;
			is_lostsight <- false;
			color <- #blue;
			//id_I<-0;
		}
		
			//Create the number of Recovered spontaneously hosts 
		create Host number: R2_people {
			is_susceptible <- false;
			is_infected <- false;
			is_latent_e <- false;
			is_latent_f <- false;
			is_recovered_1 <- false;
			is_recovered_2<- true;
			is_transfered <- false;
			is_lostsight <- false;
			color <- rgb(#77B5FE);
		}
		
			//Create the number of Transfered hosts 
		create Host number: T_people {
			is_susceptible <- false;
			is_infected <- false;
			is_latent_e <- false;
			is_latent_f <- false;
			is_recovered_1 <- false;
			is_recovered_2<- false;
			is_transfered <- true;
			is_lostsight <- false;
			color <- #gray;
		}
		
			//Create the number of Lost Sight hosts 
		create Host number: K_people {
			is_susceptible <- false;
			is_infected <- false;
			is_latent_e <- false;
			is_latent_f <- false;
			is_recovered_1 <- false;
			is_recovered_2<- false;
			is_transfered <- false;
			is_lostsight <- true;
			color <- #magenta;
		}
		
	}


	reflex stop_simulation when: (cycle = 50) {
					do pause ;
			}
}
//Grid that will be used to discretize space
grid model_grid width: 50 height: 50 {
		rgb color <- #BDD9FD;
		list<model_grid> neighbours <- (self neighbors_at range) of_species model_grid;
	}
	
//Species host which represents the host of the disease
species Host {
	
	//Different booleans to know in which state is the host
	bool is_susceptible <- true;
	bool is_infected <- false;
	bool is_latent_e <- false;
	bool is_latent_f <- false;
	bool is_recovered_1 <- false;
	bool is_recovered_2<- false;
	bool is_transfered <- false;
	bool is_lostsight <- false;
	
	//float is_infected_sinse<- time;
	
	//Color of the host
	rgb color <- #green;
	
	//Location of the agent among the grid
	model_grid myPlace;
	
	
	//Count of neighbors infected 
    int ngb_infected_number function: 
    self neighbors_at(range) count(each.is_infected);
    
	
	init {
		//The location is chosen randomly
		myPlace <- one_of(model_grid);
		set location <- myPlace.location;
	}
	//Reflex to move the agent in the neighbours cells
	reflex basic_move {
		myPlace <- one_of(myPlace.neighbours); 
		location <- myPlace.location;

// ************* debut pour Is_Susceptible **************	1	 
	}
//Recruter les individus

	reflex recruter_S when: is_susceptible {
		create Host number: Gamma*N{
			bool is_susceptible <- true;
			bool is_infected <- false;
			bool is_latent_e <- false;
			bool is_latent_f <- false;
			bool is_recovered_1 <- false;
			bool is_recovered_2<- false;
			bool is_transfered <- false;
			bool is_lostsight <- false;
			rgb color <- #green;
			
		} 
		
	}
	
	reflex natural_die when: ((is_susceptible or is_recovered_1 or is_recovered_2 or is_latent_e or is_latent_f) and flip(mu1)){
        do die;
    }
    
    reflex natural_and_infection_die when: ((is_infected or is_transfered or is_lostsight) and flip(mu1)) or ((is_infected or is_transfered or is_lostsight) and flip(mu2)){
        do die;
    }

			
	reflex save_data when: every(step){
		//save the following text into the given text file. Note that each time the save statement is used, a new line is added at the end of the file.
		//save (""+cycle+ ";"+first(math_TB).S+";"+first(math_TB).I+";"+first(math_TB).Le+";"+first(math_TB).Lf+";"+first(math_TB).T+";"+first(math_TB).K+";"+first(math_TB).R1+";"+first(math_TB).R2) to: "../results/dataall.txt" rewrite: false;
		//save (""+cycle+";"+first(math_TB).S+";"+first(math_TB).I+first(math_TB).Le+first(math_TB).Lf+first(math_TB).T+first(math_TB).K+";"+first(math_TB).R1+first(math_TB).R2) to: "../results/datasir.txt" rewrite: false;
		//save [name,speed, size] to: "../results/bug.csv" type:"csv"
		//save [cycle,first(math_TB).I/100] to: "../results/dataI.csv" type:"csv" rewrite: false ;
		//save the following text into the given text file. Note that each time the save statement is used, a new line is added at the end of the file.
		save [cycle, (Host as list) count (each.is_susceptible),(Host as list) count (each.is_infected),(Host as list) count (each.is_latent_e),(Host as list) count (each.is_latent_f),(Host as list) count (each.is_transfered),(Host as list) count (each.is_lostsight),(Host as list) count (each.is_recovered_1),(Host as list) count (each.is_recovered_2)] to: "../results/datas.csv" type:"csv" rewrite: false ;
		//save [cycle, nb_susceptible, nb_infected, nb_latent_e, nb_latent_f, nb_recovered_1, nb_recovered_2, nb_transfered, nb_lostsight] to: "../results/datasirdfe.csv" type:"csv" rewrite: false;
		
		
		
	}
	//Reflex to pass the agent to the state infected 
	reflex become_infected_from_S when: is_susceptible {
			//Probability of being infected according to the number of infected among the neighbours 
    		if (flip(1 - (1 - lambda)  ^ ngb_infected_number))  {  //alpha ==> ngb_infected_number
        		is_susceptible <- false;
				is_infected <- true;
				is_latent_e <- false;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #red; 
				  			
			}    				
	}
	
	//Reflex to pass the agent to the state Latent Early
		reflex become_latent_e_from_S when: is_susceptible{ 
			if (flip(1 - (1 - lambda)  ^ ngb_infected_number))  {
		is_susceptible <- false;
		is_infected <- false;
		is_latent_e <- true;
		is_latent_f <- false;
		is_recovered_1 <- false;
		is_recovered_2<- false;
		is_transfered <- false;
		is_lostsight <- false;
		color <- #orange;
	} 

}




// ************* fin pour Is_Susceptibles **************	


// ************* debut pour Is_latent_e **************	2		
	
	//Reflex to pass the agent to the state recovered_1
	reflex become_recovered_2_from_Le when: (is_latent_e  and flip(g1)) {
		is_susceptible <- false;
		is_infected <- false;
		is_latent_e <- false;
		is_latent_f <- false;
		is_recovered_1 <- true;
		is_recovered_2<- false;
		is_transfered <- false;
		is_lostsight <- false;
		color <- rgb(#77B5FE);
	} 
	
	//Reflex to pass the agent to the state recovered_1
	reflex become_recovered_1_from_Le when: (is_latent_e  and flip(g2)) {
		is_susceptible <- false;
		is_infected <- false;
		is_latent_e <- false;
		is_latent_f <- false;
		is_recovered_1 <- false;
		is_recovered_2<- true;
		is_transfered <- false;
		is_lostsight <- false;
		color <- #blue;
	} 
	
		//Reflex to pass the agent to the state infected
	reflex become_infected_from_Le when: (is_latent_e  and flip(q)) {
		is_susceptible <- false;
		is_infected <- true;
		is_latent_e <- false;
		is_latent_f <- false;
		is_recovered_1 <- false;
		is_recovered_2<- false;
		is_transfered <- false;
		is_lostsight <- false;
		color <- #red;
	} 

// ************* fin pour Is_latent_e **************


// ************* debut pour Is_latent_f ************** 3

	reflex become_recovered_2_from_Lf when: (is_latent_f  and flip(k1)) {
		is_susceptible <- false;
		is_infected <- false;
		is_latent_e <- false;
		is_latent_f <- false;
		is_recovered_1 <- true;
		is_recovered_2<- false;
		is_transfered <- false;
		is_lostsight <- false;
		color <- #blue;
	} 
	
		reflex become_recovered_1_from_Lf when: (is_latent_f  and flip(k2)) {
		is_susceptible <- false;
		is_infected <- false;
		is_latent_e <- false;
		is_latent_f <- false;
		is_recovered_1 <- false;
		is_recovered_2<- true;
		is_transfered <- false;
		is_lostsight <- false;
		color <- rgb(#77B5FE) ;
	} 
	
		reflex become_infeted_from_Lf when: (is_latent_f  and flip(w)) {
		is_susceptible <- false;
		is_infected <- true;
		is_latent_e <- false;
		is_latent_f <- false;
		is_recovered_1 <- false;
		is_recovered_2<- false;
		is_transfered <- false;
		is_lostsight <- false;
		color <- #red;
	} 

// ************* fin pour Is_latent_f **************			
	
// ************* debut pour Is_Infected **************	4

	reflex become_recovered_1_from_I  when: (is_infected  and flip(gamma)) { 
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- false;
				is_latent_f <- false;
				is_recovered_1 <- true;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #blue;       			
			}  

					  				
	reflex become_transferd_from_I when: (is_infected  and flip(beta)) { 
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- false;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- true;
				is_lostsight <- false;
				color <- #gray;       			
			} 	
	
	reflex become_lostsight_from_I when: (is_infected  and flip(v)) { 
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- false;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- true;
				color <- #magenta;       			
			}
		reflex become_recovered_2_from_I  when: (is_infected  and flip(sigma)) { 
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- false;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- true;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #blue;       			
			}	
				
	reflex become_latent_e_from_I when: (is_infected  and flip(r1)) { 
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- true;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #orange;       			
			} 
			
					
	
//*************************fin pour Is_Infected 

//*************************debut pour Is_Transfered 5
		  				
	reflex become_latent_e_from_T when: (is_transfered  and flip(r2)) { 
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- true;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #orange;       			
			} 

//*************************fin pour Is_Transfered 

//*************************debut pour Is_Transfered 6
		  				
	reflex become_latent_e_from_K when: (is_lostsight and flip(r3)) { 
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- true;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #orange;       			
			} 

//*************************fin pour Is_Transfered 
		
//*************************debut pour Is_Recovered_1  --  7

reflex become_infected_from_R1 when: is_recovered_1 {
			//Probability of being infected according to the number of infected among the neighbours
    		if (flip(1- (1 - lambda) ^ ngb_infected_number)) {  //alpha ==> ngb_infected_number
        		is_susceptible <- false;
				is_infected <- true;
				is_latent_e <- false;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #red; 	
				//id_I<- id_I + 1;			  			
			}    				
	}

reflex become_latent_e_from_R1 when: is_recovered_1 {
			//Probability of being infected according to the number of infected among the neighbours
    		if (flip(1 - (1 - lambda)  ^ ngb_infected_number)) {  //alpha ==> ngb_infected_number
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- true;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #orange; 
				  			
			}    				
	}

//*************************debut pour Is_Recovered_1  --  

//*************************debut pour Is_Recovered_2  --  8

reflex become_infected_from_R2 when: is_recovered_2 {
			//Probability of being infected according to the number of infected among the neighbours
    		if (flip(1 - (1 - lambda) ^ ngb_infected_number ^ r )) {  //alpha ==> ngb_infected_number
        		is_susceptible <- false;
				is_infected <- true;
				is_latent_e <- false;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #red; 
				  			
			}    				
	}

reflex become_latent_e_from_R2 when: is_recovered_2 {
			//Probability of being infected according to the number of infected among the neighbours
    		if (flip(1 - (1 - lambda)  ^ ngb_infected_number )) {  //alpha ==> ngb_infected_number
        		is_susceptible <- false;
				is_infected <- false;
				is_latent_e <- true;
				is_latent_f <- false;
				is_recovered_1 <- false;
				is_recovered_2<- false;
				is_transfered <- false;
				is_lostsight <- false;
				color <- #orange; 
			}    				
			
	}


//*************************fin pour Is_Recovered_2  -- 

	reflex compter when: cycle >1{
		//write("Les Susceptibles   "+S_people); 
		//write("ngb_infected_number "+ngb_infected_number);  
		//write("(1 - (1 - beta) et ngb_inf  "+(flip(1 - (1 - beta)  ^ ngb_infected_number)));
		//write("(1 - (1 - beta)  "+(1 - (1 - beta)));  
		//write("FLIP   "+flip(gamma)); 
		
		/*write("TOTAL"+((Host as list) count (each.is_susceptible) + 
			(Host as list) count (each.is_infected) + 
			(Host as list) count (each.is_latent_e) + 
			(Host as list) count (each.is_latent_f) +
			(Host as list) count (each.is_transfered) +
			(Host as list) count (each.is_lostsight) +
			(Host as list) count (each.is_recovered_1) +
			(Host as list) count (each.is_recovered_2)
		) ); */
		
		//write("TOTAL DES SUSCEPTIBLES"+((Host as list) count (each.is_susceptible)) );
		int jour<-#day;
		write "temps : "+ jour;
		write "cycle : "+ cycle;
		write "ngb_infected_number : "+ ngb_infected_number;
		write "valeurzflip : "+flip(1 - (1 - lambda)  ^ ngb_infected_number);
		write "rangerange : "+range;
		
	}


	aspect basic {
		draw circle(0.6) color: color;
	}

}
experiment ABM_TB type: gui {	
	
	//Number of all hosts
	int nb_hosts <- N update: length(Host);
	//Number of Susceptible hosts
	int nb_susceptible <- S_people update:  Host count (each.is_susceptible);
	//Number of infected hosts
	int nb_infected <- I_people update:  Host count (each.is_infected);
	//Number of latent_e hosts
	int nb_latent_e <- Le_people update:  Host count (each.is_latent_e);
	//Number of latent_f hosts
	int nb_latent_f <- Lf_people update:  Host count (each.is_latent_f);
	//Number of recovered_1 hosts
	int nb_recovered_1 <- R1_people update:  Host count (each.is_recovered_1);
	//Number of recovered_2 hosts
	int nb_recovered_2<- R2_people update:  Host count (each.is_recovered_2);
	//Number of transfered hosts
	int nb_transfered <- T_people update:  Host count (each.is_transfered);
	//Number of lostsight hosts
	int nb_lostsight<- K_people update:  Host count (each.is_lostsight);

	//population
	parameter 'Number of Susceptible: S' type: int var: S_people category: "Initial population";
	parameter 'Number of Active TB: I' type: int var: I_people category: "Initial population";
	parameter 'Number early latent TB: Le' type: int var: Le_people category: "Initial population";
	parameter 'Number late latent TB: Lf' type: int var: Lf_people category: "Initial population";
	parameter 'Number of transfered TB people: T' type: int var: T_people category: "Initial population";
	parameter 'Number of TB  interrupt treatment: K' type: int var: K_people category: "Initial population";
	parameter 'Number of Recovered after treatment: R1' type: int var: R1_people category: "Initial population";
	parameter 'Number of spontaneously Recovered: R2' type: int var: R2_people category: "Initial population";
	//my parameters   
	parameter 'Range of contamination' type: int var: range<- 5 category: "Parameters"; //Range 
	parameter 'Recrutment: Gamma' type: float var: Gamma <- 0.000111 category: "Parameters"; //Recrutment
	parameter 'Naturaly mortality: mu1' type: float var: mu1 <- 0.00222 category: "Parameters";//rate of mortality related to TB infection  
	parameter 'Mortality due to TB: mu2' type: float var: mu2 <- 0.04  category: "Parameters"; //rate of mortality not related to TB infection 0.0003
	parameter 'Rate of transfer: beta (I->T)' type: float var:beta <- 0.009 category: "Parameters"; //rate of transferred people in a other hospital
	parameter 'Rate of recovered after treatment: gamma (I->R1)' type: float var:gamma <- 0.84 category: "Parameters"; //rate of recovered after treatment process
	parameter 'Rate of spontaneously recovered: sigma (I->R2)' type: float var: sigma<-0.10 category: "Parameters"; //rate of spontaneously recovered
	parameter 'Probability of transmission: lambda' type: float var: lambda <- 0.10 category: "Parameters"; // rate of transmission
	parameter 'Fast progressor to ATB: q (Le->I)' type: float var: q <-1.5  category: "Parameters"; // Rate of progression to ATB: q (Le->I) per year
	parameter 'Rate of progression to ATB: w (Lf->I)' type: float var: w <-0.0005 category: "Parameters"; // Rate of progression to ATB: w (Lf->I)
	parameter 'Slow progression to Latent late: h (Le->Lf)' type: float var:h <- 0.05  category: "Parameters"; //rate of re-infection from R1 and R2 to Le
	parameter 'Rate of re-infection: r1 (I->Le)' type: float var: r1 <- 0.0002  category: "Parameters"; //rate of re-infection from I to Le
	parameter 'Rate of re-infection: r2 (T->Le)' type: float var: r2 <- 0.002  category: "Parameters"; //rate of re-infection from T to Le
	parameter 'Rate of re-infection: r3 (K->Le)' type: float var: r3 <- 0.002  category: "Parameters"; //rate of re-infection from K to Le
	parameter 'Recovered after treatment: g1 (Le->R1)' type: float var: g1 <- 0.00  category: "Parameters"; //rate of recovered after treatment process from Le to R1
	parameter 'Recovered spontaneously: g2 (Le->R2)' type: float var: g2 <- 0.0  category: "Parameters"; // rate of spontaneously recovered from Le to R2
	parameter 'Recovered after treatment: k1 (Lf->R1)' type: float var: k1 <- 0.00 category: "Parameters"; //rate of recovered after treatment process from Lf to R1
	parameter 'Recovered spontaneously: k2(Lf->R2)' type: float var: k2 <- 0.00  category: "Parameters"; // rate of spontaneously recovered from Lf to R2
	parameter 'Rate of reactivation: r (R1 & R2 -> Le)' type: float var: r <- 0.0436  category: "Parameters"; // Rate of reactivation: r (R1 & R2 -> Le)
	parameter'Probability lost to follow up while on treatment: v (I->K)' type: float var: v <- 0.06  category: "Parameters"; //Rate of treatment interuption v (I->K)
	 
	output {
			display sir_display { 
			grid model_grid lines: #gray;
			species Host aspect: basic;	
		}
	
display TB_ALL{ 
			chart "" type: series background: #white 
			x_label: "Time (Years)" y_label: "Number of S(t), Le(t), Lf(t), I(t), K(t), R1(t), R2(t) and T(t), " 
			tick_font: 'Times New Roman' tick_font_size: 14 tick_font_style: 'plain'
			legend_font_size: 18
			label_font: 'Arial' label_font_size: 14 label_font_style: 'plain'{
				data 'S' value: (Host as list) count (each.is_susceptible) color: #green;
				data 'Le' value: (Host as list) count (each.is_latent_e) color: #orange;
				data 'Lf' value: (Host as list) count (each.is_latent_f) color: #yellow;
				data 'I' value: (Host as list) count (each.is_infected) color: #red;
				data 'K' value: (Host as list) count (each.is_lostsight) color: #magenta;
				data 'R1' value: (Host as list) count (each.is_recovered_1) color: #blue;
				data 'R2' value: (Host as list) count (each.is_recovered_2) color: rgb(#77B5FE);				
				data 'T' value: (Host as list) count (each.is_transfered) color: #gray;
			}
		}
		

	display ABM { 
			chart '' type: series background: #black style: exploded {
				data 'S' value: (Host as list) count (each.is_susceptible) color: #green;
				data 'Le' value: (Host as list) count (each.is_latent_e) color: #orange;
				data 'Lf' value: (Host as list) count (each.is_latent_f) color: #yellow;
				data 'I' value: (Host as list) count (each.is_infected) color: #red;
				data 'K' value: (Host as list) count (each.is_lostsight) color: #magenta;
				data 'R1' value: (Host as list) count (each.is_recovered_1) color: #blue;
				data 'R2' value: (Host as list) count (each.is_recovered_2) color: rgb(#77B5FE);				
				data 'T' value: (Host as list) count (each.is_transfered) color: #gray;
				
			}
		}
		
		display Camember { 
			chart '' type: pie style: "3d" {
				data 'S' value: (Host as list) count (each.is_susceptible)color: #green;
				data 'Le' value: (Host as list) count (each.is_latent_e) color: #orange;
				data 'Lf' value: (Host as list) count (each.is_latent_f) color: #yellow;
				data 'I' value: (Host as list) count (each.is_infected) color: #red;
				data 'K' value: (Host as list) count (each.is_lostsight) color: #magenta;
				data 'R1' value: (Host as list) count (each.is_recovered_1) color: #blue;
				data 'R2' value: (Host as list) count (each.is_recovered_2) color: rgb(#77B5FE);				
				data 'T' value: (Host as list) count (each.is_transfered) color: #gray;
				
			}
		}
 }
}
