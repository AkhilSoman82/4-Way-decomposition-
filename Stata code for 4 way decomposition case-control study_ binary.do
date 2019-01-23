****** Run code from the line below till the end at a single stretch *******

**Start of code **

set varabbrev off, perm
cap prog drop calc2
prog calc2, rclass

logit Y A##M C1 C2 C3 C4 // Outcome model

scalar t1=_b[1.A]
scalar t2=_b[1.M]
scalar t3=_b[1.A#1.M]

logit M A C1 C2 C3 C4 if Y==0 // Mediator model fit among among controls 

scalar b0=_b[_cons]
scalar b1=_b[A]

//Covariates
scalar bc1 = _b[C1]
scalar bc2 = _b[C2]
scalar bc3 = _b[1.C3]
scalar bc4 = _b[1.C4]

// Assigning level of covariates // for continuous covariates, just take the mean of // their distribution in full sample 
sum C1 
scalar cc1= r(mean)
sum C2 
scalar cc2= r(mean)
scalar cc3=1 // at level 1 of binary covariate C3. Any level can be assigned based on // requirement                
scalar cc4=3 // at level 3 of 3 category covariate C4

// Calculating bcc - Calculating sum of products of coefficients of covariates from 
// mediator model and level of covariate

scalar bcc = bc1*cc1 + bc2*cc2 + bc3*cc3+ bc4*cc4 
// Additional values assigned
scalar a1=1 // level 1 of exposure A
scalar a0=0 // level 0 of exposure A
scalar m0=0 // level 0 of binary mediator
scalar mstar=0 // level of mediator at which CDE is calculated



// 2-way decomposition or mediation – Calculating coefficients for natural direct effect (NDE), natural indirect effect (NIE) and total effect (total)

scalar lnde    = ln((exp(t1*a1)*(1+exp(t2+t3*a1+b0+b1*a0+bcc))) ///     
                  /(exp(t1*a0)*(1+exp(t2+t3*a0+b0+b1*a0+bcc))))
scalar lnie    = ln(((1+exp(b0+b1*a0+bcc))*(1+exp(t2+t3*a1+b0+b1*a1+bcc))) ///  
                  /((1+exp(b0+b1*a1+bcc))* (1+exp(t2+t3*a1+b0+b1*a0+bcc)))) 
scalar ltotal  = ln((exp(t1*a1)*(1+exp(b0+b1*a0+bcc))* ///
                 (1+exp(b0+b1*a1+bcc+t2+t3*a1)))/(exp(t1*a0)* /// 
                 (1+exp(b0+b1*a1+bcc))* (1+exp(b0+b1*a0+bcc+t2+t3*a0))))

// 4- way decomposition - Calculating coefficients for controlled direct effect (CDE), 
// reference interaction(INTref), mediated interaction(INTmed), pure indirect effect 
// (PIE)

scalar lcde    = ln(exp(t1 + t3*mstar)*(a1-a0))

scalar lIntref = ln((exp(t1*(a1-a0)-t2*mstar-t3*a0*mstar)* ///  
                 (1+exp(b0+b1*a0+bcc+t2+t3*a1))) /(1+exp(b0+b1*a0+bcc)) ///
                 -(exp(-t2*mstar-t3*a0*mstar)*(1+exp(b0+b1*a0+bcc+t2+t3*a0))) ///  
                 /(1+exp(b0+b1*a0+bcc)) - exp((t1+t3*mstar)*(a1-a0)) + 1)
												 
scalar lIntmed = ln((exp(t1*(a1-a0)-t2*mstar-t3*a0*mstar)* /// 
                 (1+exp(b0+b1*a1+bcc+t2+t3*a1)) /(1+exp(b0+b1*a1+bcc))) ///
                 - (exp(-t2*mstar-t3*a0*mstar)*(1+exp(b0+b1*a1+bcc+t2+t3*a0)) ///  
                 /(1+exp(b0+b1*a1+bcc)))- exp(t1*(a1-a0)-t2*mstar-t3*a0*mstar) /// 
                 *(1+exp(b0+b1*a0+bcc+t2+t3*a1))/(1+exp(b0+b1*a0+bcc)) ///
		    + exp(-t2*mstar-t3*a0*m0)*(1+exp(b0+b1*a0+bcc+t2+t3*a0)) /// 
                  /(1+exp(b0+b1*a0+bcc)))			
		
scalar lpie    = ln((1+exp(b0+b1*a0+bcc))*(1+exp(b0+b1*a1+bcc+t2+t3*a0)) ///
                 /((1 + exp(b0+b1*a1+bcc))*(1+exp(b0+b1*a0+bcc+t2+t3*a0)))) 

// Calculating coefficients for each 4-way component and total effect
scalar cde_comp    = (exp(t1*(a1-a0)+t2*mstar+t3*a1*mstar)*(1+exp(b0+b1*a0+bcc))/ ///
                     (1+exp(b0+b1*a0+bcc+t2+t3*a0)))-(exp(t2*mstar+t3*a0*mstar)* ///    
                     (1+exp(b0+b1*a0+bcc))/(1+exp(b0+b1*a0+bcc+t2+t3*a0)))		
scalar INTref_comp = exp(t1*(a1-a0))*(1+exp(b0+b1*a0+bcc+t2+t3*a1)) ///  
                     /(1+exp(b0+b1*a0+bcc+t2+t3*a0)) - (1) ///
                     -exp(t1*(a1-a0)+t2*mstar+t3*a1*mstar)*(1+exp(b0+b1*a0+bcc)) /// 
                     /(1+exp(b0+b1*a0+bcc+t2+t3*a0))+ exp(t2*mstar+t3*a0*mstar)* /// 
                     (1+exp(b0+b1*a0+bcc))/(1+exp(b0+b1*a0+bcc+t2+t3*a0))
scalar INTmed_comp = exp(t1*(a1-a0))*(1+exp(b0+b1*a1+bcc+t2+t3*a1))* /// 
                     (1+exp(b0+b1*a0+bcc))/((1+exp(b0+b1*a0+bcc+t2+t3*a0)) ///
                     *(1+exp(b0+b1*a1+bcc)))-(1+exp(b0+b1*a1+bcc+t2+t3*a0))* ///  
                     (1+exp(b0+b1*a0+bcc)) /((1+exp(b0+b1*a0+bcc+t2+t3*a0))* ///    
                     (1+exp(b0+b1*a1+bcc))) - exp(t1*(a1-a0))* ///  
                     (1+exp(b0+b1*a0+bcc+t2+t3*a1)) ///
                      /(1+exp(b0+b1*a0+bcc+t2+t3*a0)) + (1)  
scalar pie_comp    = (1+exp(b0+b1*a0+bcc))*(1+exp(b0+b1*a1+bcc+t2+t3*a0)) ///
                      /((1 + exp(b0+b1*a1+bcc))*(1+exp(b0+b1*a0+bcc+t2+t3*a0))) -(1)
scalar total       = (exp(t1*a1)*(1+exp(b0+b1*a0+bcc))* ///
                     (1+exp(b0+b1*a1+bcc+t2+t3*a1))) /(exp(t1*a0)* ///    
                     (1+exp(b0+b1*a1+bcc))* (1+exp(b0+b1*a0+bcc+t2+t3*a0)))


// Retrieving the values of each coefficient calculated above
return scalar lnie=lnie
return scalar lnde=lnde
return scalar ltotal=ltotal						
return scalar lcde=lcde	
return scalar lIntref=lIntref	
return scalar lIntmed=lIntmed
return scalar lpie=lpie		
return scalar cde_comp = cde_comp
return scalar INTref_comp = INTref_comp
return scalar INTmed_comp = INTmed_comp
return scalar pie_comp  = pie_comp 
return scalar total=total
// Calculating value for total excess relative risk (terr)
scalar terr   = cde_comp +INTref_comp + INTmed_comp + pie_comp  
// Calculating the values for each of the 4 components of the total excess risk
scalar errCDE                    = cde_comp*(total-1)/terr
scalar errINTref 			= INTref_comp*(total-1)/terr
scalar errINTmed 			= INTmed_comp*(total-1)/terr
scalar errPIE    			= pie_comp*(total-1)/terr
// Assigning the values for proportions of total excess risk that is due to each of  
// the component
scalar PropCDE    			= cde_comp/terr
scalar PropINTref 			= INTref_comp/terr
scalar PropINTmed 			= INTmed_comp/terr
scalar PropPIE    			= pie_comp/terr

// Assigning the values of overall proportions of risk attributable to mediation,  
// interaction, and proportion eliminated
scalar PropMediated 	      = (pie_comp+INTmed_comp)/terr
scalar PropAttribInteraction     = (INTref_comp+INTmed_comp)/terr
scalar PropEliminated            = (INTref_comp+INTmed_comp+pie_comp)/terr				
// Retrieving the values for each of the above
return scalar terr               = terr
return scalar errCDE   		= errCDE 
return scalar errINTref 		= errINTref
return scalar errINTmed 		= errINTmed 
return scalar errPIE    		= pie_comp
return scalar PropCDE    		= PropCDE  
return scalar PropINTref 		= PropINTref 
return scalar PropINTmed 		= PropINTmed 
return scalar PropPIE    		= PropPIE 
return scalar PropMediated 			= PropMediated
return scalar PropAttribInteraction = PropAttribInteraction
return scalar PropEliminated        = PropEliminated 						
end 
calc2
return list	
	
****End of code *** 

***Run the code from start till the line above****


// Running the above code will give an output of estimates (coefficients) only, of all // parameters 


**Code for calculating confidence intervals and exponentiated risk estimates
// Stata codes for - bootstrap - procedure to calculate the 95% CIs for estimate of each parameter. Any number of repetition (reps) can be assigned. Random seed number (seed) should be assigned.	


*Code: 
bootstrap lcde=r(lcde) lIntref=r(lIntref) lIntmed = r(lIntmed)  lpie=r(lpie) ///  
    ltotal=r(ltotal)  lnde=r(lnde) lnie=r(lnie), reps(2000) seed(438766) nodrop: calc2 
 
// Running the above command will create an output of results with estimates of
// observed coefficients, Bootstrap standard error, z, P>|z| and 95% CI in a table of  // rows and 6 columns 

// Steps to exponentiate the values of each observed coefficient and corresponding CIs // in the results table 
matrix T= r(table) // captures the real matrix returned by -bootstrap-
matrix list T

// Calculating the exponentiated results 
display "CDE=" exp(T[1,1]), "LB=" exp(T[5,1]), "UB=" exp(T[6,1])
display "INTref=" exp(T[1,2]), "LB=" exp(T[5,2]), "UB=" exp(T[6,2])
display "INTmed=" exp(T[1,3]), "LB=" exp(T[5,3]), "UB=" exp(T[6,3])
display "PIE="exp(T[1,4]), "LB=" exp(T[5,4]), "UB=" exp(T[6,4])
display "TE=" exp(T[1,5]), "LB=" exp(T[5,5]), "UB=" exp(T[6,5])
display "NDE=" exp(T[1,6]), "LB=" exp(T[5,6]), "UB=" exp(T[6,6])
display "NIE=" exp(T[1,7]), "LB=" exp(T[5,7]), "UB=" exp(T[6,7])

// Calculating the CIs using bootstrap for 4- way components, 4 components of excess // relative risks, and proportions 
bootstrap cde_comp = r(cde_comp) INTref_comp = r(INTref_comp) ///
          INTmed_comp = r(INTmed_comp) pie_comp = r(pie_comp) /// 
          terr = r(terr) errCDE = r(errCDE) errINTref = r(errINTref) ///
          errINTmed = r(errINTmed) errPIE = r(errPIE) PropCDE = r(PropCDE) ///    
          PropINTref = r(PropINTref) PropINTmed = r(PropINTmed) ///
          PropPIE = r(PropPIE) PropMediated = r(PropMediated) ///
          PropAttribInteraction = r(PropAttribInteraction) ///
          PropEliminated =r(PropEliminated), reps(2000) seed(438766) /// 
          saving(`boot_results') nodrop: calc2	
