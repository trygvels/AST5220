
  �$  [   k820309    P          16.0        ��X                                                                                                           
       time_mod.f90 TIME_MOD                                                     
                                                           
                                                           
                                                           
                                                                                                                                                                                                                                                                    
                   
                  )�mr�D<                                                         
                 
                 �k$	�?        8.3D-5#         @                                  	                 	   #YSTART 
   #X1    #X2    #EPS    #H1    #HMIN    #DERIVS    #RKQS    #OUTPUT !             
                               
                   
               &                                                     
                                      
                
                                      
                
                                      
                
                                      
                
                                      
      #         @                                        	               #X    #Y    #DYDX              
                                     
                
                                                   
              &                                                                                                       
               &                                           #         @                                        	            	   #Y    #DYDX    #X    #HTRY    #EPS    #YSCAL    #HDID    #HNEXT    #DERIVS              
                                                  
               &                                                     
                                                   
 	             &                                                     
                                    
                 
                                     
                
                                     
                
                                                   
 
             &                                                                                         
                                                     
       #         @                                        	               #X    #Y    #DYDX               
                                     
                
                                                   
              &                                                                                                        
               &                                           #         @                                   !     	               #X "   #Y #             
                                "     
                
                                #                   
              &                                           #         @    @                             $                 	   #Y %   #DYDX &   #X '   #HTRY (   #EPS )   #YSCAL *   #HDID +   #HNEXT ,   #DERIVS -             
                               %                   
 
              &                                                     
                                 &                   
              &                                                     
                                '     
                 
                                 (     
                
                                 )     
                
                                 *                   
              &                                                                                     +     
                                                 ,     
       #         @                                   -     	               #X .   #Y /   #DYDX 0             
                                .     
                
                                /                   
              &                                                                                    0                   
 	              &                                           #         @                                  1                    #X 2   #Y 3   #YP1 4   #YPN 5   #Y2 6             
                                2                   
              &                                                     
                                 3                   
              &                                                     
                                 4     
                
                                 5     
                                                6                   
               &                                                                                       7     
                
          )       -DT�!	@        3.141592653589793238462643383279502884197                                            8     
                 
                 u��&iW�=        6.67258D-11                                            9     
                 
                 y�&1��?        0.224D0                                            :     
                   
                  �<kM$ˆ:                                                    ;     
                 
                 Zd;�O��?        0.046D0                                            <     
                   
                  �ky�z[�?                                                    =     
                 
                    JxޱA        2.99792458D8                                            >     
                 
                                 0.D0%         @                               ?                    
       #XA @   #YA A   #Y2A B   #X C             
                                @                   
              &                                                     
                                 A                   
              &                                                     
                                 B                   
              &                                                     
                                 C     
                 @                                D                     @ @                              E                   
                &                                                    @                                F                   
                &                                                      @                                G                     @ @                              H                   
                &                                                    @ @                              I                   
                &                                                    @ @                              J                   
                &                                           #         @                                   K                     #         @    @                             L                    #X M   #ETA N   #DERIVATIVE O                                                           
  @                              M     
                
                                 N                   
              &                                                     D                                O                   
               &                                           #         @    @                             P                    #X Q   #Y R             
                                 Q     
                
                                 R                   
              &                                           %         @                               S                    
       #X_IN T             
  @                              T     
      %         @                               U                    
       #X V             
                                 V     
      %         @                                W                    
       #X X             
                                 X     
      %         @                                Y                    
       #X Z             
                                 Z     
         �         fn#fn    �   @   J   HEALPIX_TYPES    �   @   J   PARAMS    >  @   J   SPLINE_1D_MOD    ~  @   J   ODE_SOLVER "   �  p       I4B+HEALPIX_TYPES !   .  p       DP+HEALPIX_TYPES    �  p       H_0+PARAMS      v       OMEGA_R+PARAMS "   �  �       ODEINT+ODE_SOLVER )   %  �   a   ODEINT%YSTART+ODE_SOLVER %   �  @   a   ODEINT%X1+ODE_SOLVER %   �  @   a   ODEINT%X2+ODE_SOLVER &   1  @   a   ODEINT%EPS+ODE_SOLVER %   q  @   a   ODEINT%H1+ODE_SOLVER '   �  @   a   ODEINT%HMIN+ODE_SOLVER )   �  `      ODEINT%DERIVS+ODE_SOLVER +   Q  @   a   ODEINT%DERIVS%X+ODE_SOLVER +   �  �   a   ODEINT%DERIVS%Y+ODE_SOLVER .     �   a   ODEINT%DERIVS%DYDX+ODE_SOLVER '   �  �      ODEINT%RKQS+ODE_SOLVER )   H  �   a   ODEINT%RKQS%Y+ODE_SOLVER ,   �  �   a   ODEINT%RKQS%DYDX+ODE_SOLVER )   `	  @   a   ODEINT%RKQS%X+ODE_SOLVER ,   �	  @   a   ODEINT%RKQS%HTRY+ODE_SOLVER +   �	  @   a   ODEINT%RKQS%EPS+ODE_SOLVER -    
  �   a   ODEINT%RKQS%YSCAL+ODE_SOLVER ,   �
  @   a   ODEINT%RKQS%HDID+ODE_SOLVER -   �
  @   a   ODEINT%RKQS%HNEXT+ODE_SOLVER .   ,  `      ODEINT%RKQS%DERIVS+ODE_SOLVER 0   �  @   a   ODEINT%RKQS%DERIVS%X+ODE_SOLVER 0   �  �   a   ODEINT%RKQS%DERIVS%Y+ODE_SOLVER 3   X  �   a   ODEINT%RKQS%DERIVS%DYDX+ODE_SOLVER )   �  V      ODEINT%OUTPUT+ODE_SOLVER +   :  @   a   ODEINT%OUTPUT%X+ODE_SOLVER +   z  �   a   ODEINT%OUTPUT%Y+ODE_SOLVER      �       BSSTEP+BS_MOD     �  �   a   BSSTEP%Y+BS_MOD #   1  �   a   BSSTEP%DYDX+BS_MOD     �  @   a   BSSTEP%X+BS_MOD #   �  @   a   BSSTEP%HTRY+BS_MOD "   =  @   a   BSSTEP%EPS+BS_MOD $   }  �   a   BSSTEP%YSCAL+BS_MOD #   	  @   a   BSSTEP%HDID+BS_MOD $   I  @   a   BSSTEP%HNEXT+BS_MOD %   �  `      BSSTEP%DERIVS+BS_MOD '   �  @   a   BSSTEP%DERIVS%X+BS_MOD '   )  �   a   BSSTEP%DERIVS%Y+BS_MOD *   �  �   a   BSSTEP%DERIVS%DYDX+BS_MOD %   A  p       SPLINE+SPLINE_1D_MOD '   �  �   a   SPLINE%X+SPLINE_1D_MOD '   =  �   a   SPLINE%Y+SPLINE_1D_MOD )   �  @   a   SPLINE%YP1+SPLINE_1D_MOD )   	  @   a   SPLINE%YPN+SPLINE_1D_MOD (   I  �   a   SPLINE%Y2+SPLINE_1D_MOD !   �  �       PI+HEALPIX_TYPES    n  {       G_GRAV+PARAMS    �  w       OMEGA_M+PARAMS    `  p       RHO_C+PARAMS    �  w       OMEGA_B+PARAMS $   G  p       OMEGA_LAMBDA+PARAMS    �  |       C+PARAMS     3  t       OMEGA_NU+PARAMS %   �  p       SPLINT+SPLINE_1D_MOD (     �   a   SPLINT%XA+SPLINE_1D_MOD (   �  �   a   SPLINT%YA+SPLINE_1D_MOD )   /  �   a   SPLINT%Y2A+SPLINE_1D_MOD '   �  @   a   SPLINT%X+SPLINE_1D_MOD    �  @       N_T    ;  �       X_T    �  �       A_T    S  @       N_ETA    �  �       X_ETA      �       ETA    �  �       ETA2 $   7  H       INITIALIZE_TIME_MOD      �       DERIVS       @   a   DERIVS%X    U   �   a   DERIVS%ETA "   �   �   a   DERIVS%DERIVATIVE    m!  V       OUTPUT    �!  @   a   OUTPUT%X    "  �   a   OUTPUT%Y    �"  Z       GET_ETA    �"  @   a   GET_ETA%X_IN    )#  W       GET_H    �#  @   a   GET_H%X    �#  W       GET_H_P    $  @   a   GET_H_P%X    W$  W       GET_DH_P    �$  @   a   GET_DH_P%X 