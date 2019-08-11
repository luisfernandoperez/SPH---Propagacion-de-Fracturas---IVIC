public class Kernel_y_Derivadas{
    
    double[] kernel;
    double[][] gradiente;
    double[] laplaciano;
    
	Utilidades 		Util 		= new Utilidades();
    Parametros 		Param 		= new Parametros();
    Otras_Variables Otras_Var 	= new Otras_Variables();
    

    public Kernel_y_Derivadas
    (
		int 		i, 
		double[] 	h, 
		int[] 		vecinos, 
		double[][] 	posicion
    )
    {
        int  	dimension 	= posicion[0].length;
		int 	N 			= posicion.length;
        int 	k			= 0;
        int		p_v			= 0;
        int		largo		= vecinos.length;
        int		vecino;
        double 	d2;
        double	dx; 
        double	dy;
        double	dx2;
        double	dy2;
        double	raiz;
        double	q = 0;
        double 	sigma;
        double  h_prom;
        double k_kernel = K_Kernel();
        
		if(largo == 0) largo = 1;	//Solo por si el metodo de vecinos falla
		
		
        kernel 	   = new double[largo];
		gradiente  = new double[largo][dimension];
		laplaciano = new double[largo];
		
		
        for(int j =0; j<largo; j++)
        {
			
            vecino 			= vecinos[j];
            h_prom 			= (h[i]+h[vecino])*0.5;
            dx 				= (posicion[i][0] - posicion[vecino][0]);
            dy 				= (posicion[i][1] - posicion[vecino][1]);
            dx2 			= dx*dx;
            dy2 			= dy*dy;
            d2 				= dx2 + dy2;
            raiz 			= Math.sqrt(d2);
            q 				= raiz/h_prom;
            kernel[j]     	= Kernel(i, j, q, h_prom, dimension); 
            gradiente[j]  	= Gradiente(i, j, dimension, q, h_prom, dx, dy);
            //if(i==10 && vecinos[j] ==1) System.out.println( q  + " " +  kernel[j]);
        }
		
		
    }

	public double K_Kernel()
	{
		int  kernels = Param.kernels;
		double tolerancia = 1.0;
		if( kernels == 0)	return 1.*tolerancia;
		if( kernels == 3)	return 2.*tolerancia;
		if( kernels == 4)	return 1.*tolerancia;
		else return 1.;
	}
	
    
    // Seccion de kernels

    public double Kernel(int i, int j, double q, double h_prom, int dimension)
    {
        int  	kernels = Param.kernels;
        double  kernel	= 0;
        double	sigma 	= 0;
		double k_kernel = K_Kernel();
		double 	PI 		= Param.PI;

		if( kernels == 0)	//Gausiano
		{
			sigma = 1./(PI*h_prom*h_prom);
			
		    if(q <= k_kernel)
				kernel = sigma*Math.exp(-q*q);
				
		    else if(q >= k_kernel)
				kernel = 0;
		}

		if( kernels == 3)	//Spline Cubico
		{
				   sigma 		= 15./7.0/PI/h_prom/h_prom;
            double DOS_TERCIOS 	= 0.66666;
            double UN_SEXTO    	= 0.16666;
			double k_2 			= k_kernel*0.5;
			
		    if(q < 1.0)
				kernel = sigma*(DOS_TERCIOS - q*q + 0.5*q*q*q);
				
		    if(q >= 1.0 && q < 2.0) 
				kernel = UN_SEXTO*sigma*(2.0 -q)*(2.0 -q)*(2.0 -q);
				
		    else if(q >= 2.0) 
		    	kernel = 0.0;
        }

		if( kernels == 4)	//Wendland C4
		{
					sigma = 9./(PI*h_prom*h_prom);
			double 	q2 	  = (1-q)*(1-q);
			
		    if(q < k_kernel) 
				kernel = sigma*q2*q2*q2*(1+6*q+35/3.*q*q);
				
		    else if(q >= k_kernel)
				kernel = 0;
        }



        return kernel;
    }

	public double[] Gradiente(int i, int j, int dimension, double q, double h_prom, double dx, double dy)
	{
		double h2 = h_prom*h_prom;
		double[] gradiente = new double[dimension];

		gradiente[0] = Derivada_W_Q(q, h_prom)*Derivada_Q_X(dx, q, h_prom, dimension);
		gradiente[1] = Derivada_W_Q(q, h_prom)*Derivada_Q_X(dy, q, h_prom, dimension);
		
		return gradiente;
	}



    public double Derivada_W_Q(double q, double h_prom)
    {
        int  kernels = Param.kernels;
		double sigma = 0;
		double k_kernel = K_Kernel();
		double primera_derivada = 0;

		if( kernels == 0)	//Gausiano
		{
			sigma = 2./(Math.PI*h_prom*h_prom*h_prom);
		    if(q <= k_kernel){ 	    primera_derivada = -sigma*q*Math.exp(-q*q);}
		    else if(q >= k_kernel){ primera_derivada = 0;}
		}

		if( kernels == 3)	//Spline Cubico
		{
			sigma = 15./(7.*Math.PI*h_prom*h_prom*h_prom);
            double TRES_MEDIOS = 1.5;
			double k_2 = k_kernel*0.5;
		    if(q <= k_2) 				primera_derivada = sigma*q*(TRES_MEDIOS*q - 2);
		    if(q > k_2 && q < k_kernel) primera_derivada = -0.5*sigma*(2-q)*(2-q);
		    else if(q >= k_kernel)		primera_derivada = 0;
        }

		if( kernels == 4)	//Wendland C4
		{
		    sigma = 168./(Math.PI*h_prom*h_prom*h_prom);
			double q2 = (1-q)*(1-q);
		    	 if(q <  k_kernel){ primera_derivada = -sigma*q*q2*q2*(1-q)*(1+5*q);}
		    else if(q >= k_kernel){ primera_derivada = 0;}
        }
		return primera_derivada;
    }


    public double Derivada_Q_X(double dx, double q, double h_prom, int dimension)
    {
		if(dimension == 1) return 1/h_prom;
		else
		{
			if(q==0 && dx == 0) return 0;
			else return dx/(q*h_prom);
		} 
    }

    

    
    // Seccion Kernel adaptativo

    public double[] H(double[][] posicion, double masas, double[] densidad, double[] h_inicial, int[][] vecinos)
    {
		int N = posicion.length;
		int metodo_adaptativo = Param.metodo_adaptativo;
        double sum_log  = 0;
        double factor_g = 0;
		double lambda 	= 0;
        double k_kernel = Param.k_kernel_adaptativo;
        double epsilon  = Param.epsilon_kernel_adaptativo;
		double[] densidad_kernel = new double[N];
		double densidad_inicial = Param.densidad;
		double dimension = posicion[0].length;
        double[] h = new double[N];
		double[] kernel = new double[1];
		
		
		if(metodo_adaptativo == 0)
        {
            h = Util.Copiar_Vector(h_inicial);
        }
        
		
	    if( metodo_adaptativo == 1)
	    {

            for(int i =0; i<N; i++)
            {
				kernel = new Kernel_y_Derivadas(i, h_inicial, vecinos[i],  posicion).kernel;
				densidad_kernel[i] = Otras_Var.Densidad_Kernel(masas, kernel, vecinos[i]);
                
                lambda = Math.pow(densidad_inicial/densidad_kernel[i], dimension);          
                h[i]   = lambda*h_inicial[i];
               // if(i==0)System.out.println(vecinos[i].length + " " + densidad_kernel[i] + " " +  densidad_inicial + " " +  lambda + " " + h_inicial[i]);
            }
	    }


	    if(metodo_adaptativo == 2)
	    {
			
			for(int i =0; i<N; i++)
            {
				kernel = new Kernel_y_Derivadas(i, h_inicial, vecinos[i],  posicion).kernel;
				densidad_kernel[i] = Otras_Var.Densidad_Kernel(masas, kernel, vecinos[i]);
			}
			
            for(int i =0; i<N; i++)
            {
				sum_log  = 0;
				for(int j:vecinos[i])
				{
		            sum_log += Math.log(densidad_kernel[j]);
		        }
				factor_g = Math.exp(sum_log/vecinos[i].length);
				lambda = k_kernel*Math.pow(factor_g/densidad_kernel[i], epsilon);

                h[i]   = lambda*h_inicial[i];
            }
           // System.out.println(h[20] + " " + densidad[20]);
        }
        

        
        return h;
    }
    
	//Seccion consistencia orden 0



    public double[] Cons_Kernel(int i, int[] vecinos, double masas, double[] densidad, double[] kernel, double[][] posicion)
    {
        int largo = vecinos.length;
        int vecino = 0;
        int consistencia = Param.cons_kern;
        double[] kernel_corregido = new double[largo];
		
		

        if(consistencia == 0)
        {
            double cons_0_kernel = 0;

            for(int j =0; j<largo; j++)
            {
                vecino 		   = vecinos[j];
                cons_0_kernel += kernel[j]*masas/densidad[vecino];
            }
            

            for(int j =0; j<largo; j++)
            {
                kernel_corregido[j] = kernel[j]/cons_0_kernel;
            }
            
             //System.out.println(masas + " " + densidad[i] + " " + masas/densidad[i] + " " + cons_0_kernel + " " + kernel_corregido[0]);
            
            /* MAÑANA ACOMODAR ESTO, TIENE QUE DAR COMO VALOR MAXIMO KERNEL_CORREGIDO EL DEL KERNEL SIN CORREGIR PARA REESCALAR */
            
           // Util.Producto_Vector_Escalar(maximo_kernel/maximo_nuevo, kernel_corregido);
        }

        if(consistencia == 1)
        {
            int dimension_mas_uno = Param.N + 1;
            double[] dx = new double[largo];
            double[] dy = new double[largo];
            double cons_1_kernel = 1;
            double[][] vector_correccion = new double[largo][dimension_mas_uno];
            double[][] vector_temporal_A = new double[largo][dimension_mas_uno];
            double[][] matriz_correccion = new double[dimension_mas_uno][dimension_mas_uno];
            double[][] matriz_temporal_A = new double[dimension_mas_uno][dimension_mas_uno];

            for(int j =0; j<largo; j++)
            {
                vecino 	= vecinos[j];
                dx[j] = posicion[i][0] - posicion[vecino][0];
                dy[j] = posicion[i][1] - posicion[vecino][1];
                matriz_correccion[0][0] += kernel[j]*masas/densidad[vecino];
                matriz_correccion[0][1] += kernel[j]*masas/densidad[vecino]*dx[j];
                matriz_correccion[0][2] += kernel[j]*masas/densidad[vecino]*dy[j];
                matriz_correccion[1][0] += kernel[j]*masas/densidad[vecino]*dx[j];
                matriz_correccion[1][1] += kernel[j]*masas/densidad[vecino]*dx[j]*dx[j];
                matriz_correccion[1][2] += kernel[j]*masas/densidad[vecino]*dy[j]*dx[j];
                matriz_correccion[2][0] += kernel[j]*masas/densidad[vecino]*dy[j];
                matriz_correccion[2][1] += kernel[j]*masas/densidad[vecino]*dx[j]*dy[j];
                matriz_correccion[2][2] += kernel[j]*masas/densidad[vecino]*dy[j]*dy[j];
            }


            matriz_temporal_A = Util.Inversa_Matriz_2x2_3x3(i,matriz_correccion);

            for(int j =0; j<largo; j++)
            {
                vector_correccion[j] = Util.Producto_Matriz_Vector(matriz_temporal_A, Util.x_unitario(dimension_mas_uno));
            }

            for(int j =0; j<largo; j++)
            {
                cons_1_kernel = vector_correccion[j][0] + vector_correccion[j][1]*dx[j] + vector_correccion[j][2]*dy[j];
                if(cons_1_kernel != cons_1_kernel || cons_1_kernel < 1) cons_1_kernel = 1;
                kernel_corregido[j] = cons_1_kernel*kernel[j];
            }
            //if(i==0)System.out.println(funcion_corregida + " " +  funcion + " " + matriz_temporal_A[0][0]);
        }

        if(consistencia != 0 && consistencia != 1)
        {
            for(int j =0; j<largo; j++)
            {
                kernel_corregido[j] = kernel[j];
            }
        }


        return kernel_corregido;
    }

	public double[][] Cons_Grad(int i, int[] vecinos, double masas, double[] densidad, double[][] posiciones, double[] kernel, double[] kernel_cons, double[][] grad_kernel)
	{
        int vecino = 0;
		int largo = vecinos.length;
		int dimension = posiciones[0].length;
		double tolerancia = 0.0001;
		boolean consistencia_orden_0 = Param.cons_grad;
		double[] cons_0_kernel = new double[largo];
        double[][] cons_0_grad 	   = new double[dimension][dimension];
        double[][] cons_0_grad_inv = new double[dimension][dimension];
		double[][] gradiente_sheppard = new double[largo][dimension];

		if(consistencia_orden_0)
		{
			for(int j =0; j<largo; j++)
			{
				vecino 				   = vecinos[j];
				cons_0_grad_inv[0][0] += (posiciones[vecino][0]-posiciones[i][0])*grad_kernel[j][0]*masas/densidad[vecino];
				cons_0_grad_inv[0][1] += (posiciones[vecino][0]-posiciones[i][0])*grad_kernel[j][1]*masas/densidad[vecino];
				cons_0_grad_inv[1][0] += (posiciones[vecino][1]-posiciones[i][1])*grad_kernel[j][0]*masas/densidad[vecino];
				cons_0_grad_inv[1][1] += (posiciones[vecino][1]-posiciones[i][1])*grad_kernel[j][1]*masas/densidad[vecino];
				
				//if(cons_0_grad_inv[0][0] != cons_0_grad_inv[0][0])System.out.println(i + " " +  vecino + " " + (posiciones[vecino][0]-posiciones[i][0])  + " " + grad_kernel[j][0] + " " + densidad[vecino]);
				
				if(largo == 1)
				{
					cons_0_grad_inv[0][0] = 1;
					cons_0_grad_inv[0][1] = 0;
					cons_0_grad_inv[1][0] = 0;
					cons_0_grad_inv[1][1] = 1;
				}
				
			}
			cons_0_grad = Util.Inversa_Matriz_2x2_3x3(i, cons_0_grad_inv);
			
			//int par = (int)(0.5*Param.nx*(Param.nx - 1));
			//if(i== 0) Util.Imprimir_Matriz(posiciones);
			
			for(int j =0; j<largo; j++)
			{
				gradiente_sheppard[j] = Util.Producto_Matriz_Vector(cons_0_grad, grad_kernel[j]); 
			}
		
		}
   
        if(!consistencia_orden_0)
        {
			for(int j =0; j<largo; j++)
			{
				gradiente_sheppard[j] = grad_kernel[j]; 
			}
		}
		return gradiente_sheppard;
	}

	//Seccion consistencia orden 1


	public Kernel_y_Derivadas(){}
}
