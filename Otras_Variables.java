public class Otras_Variables
{
	Utilidades Util = new Utilidades();
    Parametros Param = new Parametros();



	public double Densidad_Kernel(double masas, double[] kernel, int[] vecinos )
	{
		int vecino;
		double densidad = 0;
		double cons_arbitraria = 1/1.2350847057400323;
		
		for(int j = 0; j< vecinos.length; j++)
	    {
			densidad += masas*kernel[j];
		}
		
		return densidad;
	}
	
	public double Suma_Kernel(double[] kernel)
	{
		double Resultado = 0;

		for(int j = 0; j< kernel.length; j++)
	    {
			Resultado += kernel[j];
		}
		
		return Resultado;
	}

    public double Presion(double densidad)
    {
		double densidad_inicial = Param.densidad;
		double vel_sonido = Param.vel_sonido;

		return vel_sonido*vel_sonido*(densidad - densidad_inicial);
    }
    
    public double[][] Tensor_Estres(double presion, double[][] estres_deviatorico)
	{
		int dimension = estres_deviatorico.length;
		double[][] tensor_estres 	 = new double[dimension][dimension];
		double[][] matriz_temporal_A = new double[dimension][dimension];
		matriz_temporal_A = Util.Producto_Matriz_Escalar(-presion, Util.Matriz_Identidad_Rango_2(dimension));
		tensor_estres = Util.Suma_Matrices(matriz_temporal_A, estres_deviatorico);
		return tensor_estres;
	}

     public double[][] Tensor_Estres(double presion)
	{
		int dimension = Param.dimension;
		return Util.Producto_Matriz_Escalar(-presion, Util.Matriz_Identidad_Rango_2(dimension));
	}

	public double[][] Deformacion_Por_Estres_Plano(int i, int borde_de_particula, int[] vecinos, double masas, double[] densidad, double[][] elongacion, double[][] grad_kernel, boolean[]  vecinos_hook)
    {	 
		double[][] grad_posiciones	  = Grad_Vector_Hook(i, vecinos, elongacion, masas, densidad, grad_kernel, vecinos_hook);
		return Parte_Simetrica_Tensor_R2(grad_posiciones); 
	}	
		   
	public double[][] Grad_Vector_Hook(int i, int[] vecinos, double[][] vector, double masas, double[] densidad, double[][] grad_kernel, boolean[]  vecinos_hook)
	{
		int dimension = vector[0].length;
		int vecino;
		double[][] grad_vector = new double[dimension][dimension];
		
		for(int alfa = 0; alfa< dimension; alfa++)
        {
			for(int beta = 0; beta< dimension; beta++)
		    {
				for(int j = 0; j< vecinos.length; j++)
	    		{	
					vecino = vecinos[j];
					if(vecinos_hook[j])
						grad_vector[alfa][beta] += masas/densidad[vecino]*(vector[vecino][alfa]-vector[i][alfa])*grad_kernel[j][beta];
					//if(i==Param.N-1)System.out.println(grad_vector[alfa][beta] + " " + vector[vecino][alfa] + " " + vector[i][alfa]);
				}
			}
		}
		//if(i==0)Util.Imprimir_Matriz(grad_kernel);
		return grad_vector;
	}
	
	public double[][] Grad_Vector(int i, int[] vecinos, double[][] vector, double masas, double[] densidad, double[][] grad_kernel)
	{
		int dimension = vector[0].length;
		int vecino;
		double[][] grad_vector = new double[dimension][dimension];
		
		for(int alfa = 0; alfa< dimension; alfa++)
        {
			for(int beta = 0; beta< dimension; beta++)
		    {
				for(int j = 0; j< vecinos.length; j++)
	    		{	
					vecino = vecinos[j];
					grad_vector[alfa][beta] += masas/densidad[vecino]*(vector[vecino][alfa]-vector[i][alfa])*grad_kernel[j][beta];
					//if(i==0)System.out.println(grad_vector[alfa][beta] + " " + vector[vecino][alfa] + " " + vector[i][alfa]);
				}
			}
		}
		//if(i==0)Util.Imprimir_Matriz(grad_kernel);
		return grad_vector;
	}	
	
	public double Energia_Libre_Helmontz(double[][] tensor_deformacion, double[][] tensor_estres)
	{
		return 0.5*Util.Doble_Producto_Punto_Matrices(tensor_estres, tensor_deformacion);
	}
	
	

	public double[][] Parte_Simetrica_Tensor_R2(double[][] tensor)
	{
		int dimension = tensor.length;
		double[][] matriz_temporal_A  = new double[dimension][dimension];
		//double[][] matriz_temporal_B  = new double[dimension][dimension];
		
		for(int i =0; i < dimension; i++)
			for(int j =0; j < dimension; j++)
				matriz_temporal_A[i][j] = 0.5*(tensor[i][j] + tensor[j][i]);
				
		//matriz_temporal_A  = Util.Transpuesta(tensor);
		//matriz_temporal_A  = Util.Suma_Matrices(tensor, matriz_temporal_A);
		//return Util.Producto_Matriz_Escalar(0.5, matriz_temporal_A);
		return matriz_temporal_A;
	}
	
	public double[][] Parte_Antisimetrica_Tensor_R2(double[][] tensor)
	{
		int dimension = tensor.length;
		double[][] matriz_temporal  = new double[dimension][dimension];
		
		for(int i =0; i < dimension; i++)
			for(int j =0; j < dimension; j++)
				matriz_temporal[i][j] = 0.5*(tensor[i][j] - tensor[j][i]);
				
		return matriz_temporal;
	}

	public double[][]  Tensor_Deformacion_A_Partir_Estres(int i, double[][] tensor_estres)
	{
		int dimension = tensor_estres.length;
		double modulo_compresibilidad = Param.modulo_compresibilidad;
		double modulo_rigidez = Param.modulo_rigidez;
		double traza = Util.Traza_R2(tensor_estres);
		double inv_dim = 1/(double)(dimension + 1);
		double UNO_SOBRE_COMPRESIBILIDAD = 1/(9.0*modulo_compresibilidad);
		double UNO_SOBRE_RIGIDEZ = 1/(2.0*modulo_rigidez);
		double[][] tensor_deformacion = new double[dimension][dimension];
		double[][] matriz_temporal_A  = new double[dimension][dimension];
		double[][] matriz_temporal_B  = new double[dimension][dimension];

		matriz_temporal_A = Util.Producto_Matriz_Escalar(-inv_dim*traza, Util.Matriz_Identidad_Rango_2(dimension));
		matriz_temporal_B = Util.Suma_Matrices(tensor_estres, matriz_temporal_A);
		matriz_temporal_B = Util.Producto_Matriz_Escalar(UNO_SOBRE_RIGIDEZ, matriz_temporal_B);
		matriz_temporal_A = Util.Producto_Matriz_Escalar(UNO_SOBRE_COMPRESIBILIDAD*traza, Util.Matriz_Identidad_Rango_2(dimension));
		tensor_deformacion = Util.Suma_Matrices(matriz_temporal_A, matriz_temporal_B);
		return tensor_deformacion;
	}

	public double[][] Tensor_Estres_A_Partir_Deformacion(int i, double[][] tensor_deformacion)	
	{
		int dimension = tensor_deformacion.length;
		double modulo_compresibilidad = Param.modulo_compresibilidad;
		double modulo_rigidez = Param.modulo_rigidez;
		double traza = Util.Traza_R2(tensor_deformacion);
		double inv_dim = 1/(double)(dimension + 1);
		double constante_traza = (6.0*modulo_rigidez)/(3.0*modulo_compresibilidad + 4.0*modulo_rigidez)*(modulo_compresibilidad - 2*inv_dim*modulo_rigidez);
		//double constante_traza = Param.modulo_young;
		double DOS_RIGIDEZ = (2.0*modulo_rigidez);
		double[][] tensor_estres = new double[dimension][dimension];
		double[][] matriz_temporal_A  = new double[dimension][dimension];
		double[][] matriz_temporal_B  = new double[dimension][dimension];
		
		matriz_temporal_A = Util.Producto_Matriz_Escalar(traza*constante_traza, Util.Matriz_Identidad_Rango_2(dimension));
		matriz_temporal_B = Util.Producto_Matriz_Escalar(DOS_RIGIDEZ, tensor_deformacion);
		tensor_estres = Util.Suma_Matrices(matriz_temporal_A, matriz_temporal_B);
		return tensor_estres;	
	}


	public double[][] Tensor_Estres_A_Partir_Deformacion2(int i, double[][] tensor_deformacion)	
	{
		int dimension = tensor_deformacion.length;
		double DOS_RIGIDEZ 		= Param.DOS_RIGIDEZ;
		double constante_traza 	= Param.constante_traza_estres_plano;
		double traza 			= constante_traza*Util.Traza_R2(tensor_deformacion);
		
		double[][] tensor_estres = new double[dimension][dimension];
		
		for(int alfa =0; alfa < dimension; alfa++)
			for(int beta =0; beta < dimension; beta++)
			tensor_estres[alfa][beta] = traza*Util.Delta_Kronecker(alfa, beta) + DOS_RIGIDEZ*tensor_deformacion[alfa][beta];
			

		return tensor_estres;	
	}


	public double[][] Tensor_Deviatorico (int i, double[][] tensor_estres)
	{
		int dimension = tensor_estres.length;
		double inv_dim = 1/(double)(dimension + 1);
		double modulo_compresibilidad = Param.modulo_compresibilidad;
		double modulo_rigidez = Param.modulo_rigidez;
		double constante_traza = 6.0*modulo_rigidez/(3.0*modulo_compresibilidad + 4.0*modulo_rigidez);
		double traza = inv_dim*constante_traza*Util.Traza_R2(tensor_estres);
		double DOS_mu = 2.0*modulo_rigidez;
		double[][] deviatorico = new double[dimension][dimension];
		double[][] matriz_temporal_A = new double[dimension][dimension];

		matriz_temporal_A = Util.Producto_Matriz_Escalar(traza, Util.Matriz_Identidad_Rango_2(dimension));
		matriz_temporal_A = Util.Resta_Matrices(tensor_estres, matriz_temporal_A);
		deviatorico = Util.Producto_Matriz_Escalar(DOS_mu, matriz_temporal_A);
		return deviatorico;
	}



	
	public double[] Estres_2D_Ejes_Principales(double angulo, double[][] tensor_estres)
	{
        int dimension  = tensor_estres.length;
        double cos = 0, sen = 0, c2 = 0, sc = 0, s2 = 0;
        double[] ejes_princ_ten_estres = new double[dimension];

		 cos = Math.cos(angulo);
		 sen = Math.sin(angulo);
		 c2 = cos*cos;
		 sc = sen*cos;
		 s2 = sen*sen;
		 ejes_princ_ten_estres[0] = c2*tensor_estres[0][0] + sc*tensor_estres[1][0] + sc*tensor_estres[0][1] + c2*tensor_estres[1][1];
		 ejes_princ_ten_estres[1] = s2*tensor_estres[0][0] + sc*tensor_estres[1][0] + sc*tensor_estres[0][1] + s2*tensor_estres[1][1];
		 
		 return ejes_princ_ten_estres;    
	}

	public double[][] Estres_2D_Ejes_Normales(double angulo, double[] ejes_princ_estres)
	{
        int dimension  = ejes_princ_estres.length;
        double cos = 0, sen = 0, c2 = 0, sc = 0, s2 = 0;
		double[][] tensor_estres = new double[dimension][dimension];
		cos = Math.cos(angulo);
		sen = Math.sin(angulo);
		c2 = cos*cos;
		sc = sen*cos;
		s2 = sen*sen;
		tensor_estres[0][0] = c2*ejes_princ_estres[0]  + s2*ejes_princ_estres[1];
		tensor_estres[0][1] = sc*(ejes_princ_estres[0] -   ejes_princ_estres[1]);
		tensor_estres[1][0] = sc*(ejes_princ_estres[1] -   ejes_princ_estres[0]);
		tensor_estres[1][1] = s2*ejes_princ_estres[0]  + c2*ejes_princ_estres[1];
		 
		return tensor_estres;       
	}




    // Fuerza externa, por ahora constante
    public double[] Fuerza_Externa(int dimension)
    {
        double[] fuerza_externa = new double[dimension];
		double PI_SOBRE_180 = 0.017453292;
        double angulo_gravedad = Param.angulo_gravedad*PI_SOBRE_180;
        double mod_g = Math.abs(Param.g);

        fuerza_externa[0] = mod_g*Math.sin(angulo_gravedad);
        fuerza_externa[1] = mod_g*Math.cos(angulo_gravedad);

        return fuerza_externa;
    }

	public double Interpolacion_Escalar(int[] vecinos, double[] funcion, double masas, double[] densidad, double[] kernel)
	{
		int vecino;
		int largo = vecinos.length;
		double interpolacion_escalar = 0;
		//System.out.println(largo + " " + kernel.length);
		for(int j=0; j< largo; j++)
		{
			vecino = vecinos[j];
			interpolacion_escalar += masas/densidad[vecino]*funcion[vecino]*kernel[j];
		}
		return interpolacion_escalar;
	}
	
	
	public double[] Interpolacion_Vector(int[] vecinos, double[][] funcion, double masas, double[] densidad, double[] kernel)
	{
		int vecino;
		int largo = vecinos.length;
		int dimension = funcion[0].length;
		double[] interpolacion_vector = new double[dimension];
		for(int k=0; k< dimension; k++)
		{
			for(int j=0; j< largo; j++)
			{
				vecino = vecinos[j];
				interpolacion_vector[k] += masas/densidad[vecino]*funcion[vecino][k]*kernel[j];
			}
		}
		return interpolacion_vector;
	}

	public double[] Promedio_Vector_Vecinos(int[] vecinos, double[][] funcion)
	{
		int vecino;
		int largo = vecinos.length;
		int dimension = funcion[0].length;
		double[] interpolacion_vector = new double[dimension];
		for(int k=0; k< dimension; k++)
		{
			for(int j=0; j< largo; j++)
			{
				vecino = vecinos[j];
				interpolacion_vector[k] += funcion[vecino][k];
			}
			interpolacion_vector[k] /= largo;
		}
		return interpolacion_vector;
	}
	
	public double Promedio_Escalar_Vecinos(int[] vecinos, double[] funcion)
	{
		int vecino;
		int largo = vecinos.length;
		double interpolacion_escalar = 0;

		for(int j=0; j< largo; j++)
		{
			vecino = vecinos[j];
			interpolacion_escalar += funcion[vecino];
		}
		
		return interpolacion_escalar/largo;
	}
	
	
	public double[] Vector_Centro_de_Masa(double[][] vector, double[] masas)
	{
		int N = vector.length;
		int dimension = vector[0].length;
		double masa_total =0;
		double[] suma_parcial = new double[dimension];
		double[] resultado    = new double[dimension];

		for(int j =0; j < dimension; j++)
		{		
			for(int i =0; i < N; i++)
			{
				suma_parcial[j] += vector[i][j]*masas[i];	
				
				if(j==0)
				{
					masa_total += masas[i];
				}
			}

			resultado[j] = suma_parcial[j]/masa_total;				
		}
		return resultado;		
	}

	public double[] Elongaciones_Deformacion_Constante(double separacion , double[] tamano, double[] posicion_sin_deformar, double[][] tensor_deformacion)
	{
		int dimension = posicion_sin_deformar.length;
		double[] elongacion = new double[dimension];

		if(dimension == 1)
			elongacion[0] = tensor_deformacion[0][0]*posicion_sin_deformar[0];

		if(dimension == 2)
		{
			//elongacion[0] = tensor_deformacion[0][0]*(posicion_sin_deformar[0] + 0.5*tamano[0] ) + (2*tensor_deformacion[0][1])*(posicion_sin_deformar[1]  + tamano[1] );
			//elongacion[1] = tensor_deformacion[1][1]*(posicion_sin_deformar[1] + 0.5*tamano[1] ) + (2*tensor_deformacion[1][0])*(posicion_sin_deformar[0]  + tamano[0] );
			elongacion[0] = tensor_deformacion[0][0]*(posicion_sin_deformar[0] + 0.5*(tamano[0] - separacion)) + (2*tensor_deformacion[0][1])*(posicion_sin_deformar[1]  + 0.5*tamano[1] );
			elongacion[1] = tensor_deformacion[1][1]*(posicion_sin_deformar[1] + 0.5*(tamano[1] - separacion)) + (2*tensor_deformacion[1][0])*(posicion_sin_deformar[0]  + 0.5*tamano[0] );
			//if(i==0) System.out.println(elongacion[0] + " " + tensor_deformacion[0][0] + " " + tensor_deformacion[0][1] + " " + posicion_sin_deformar[0]);
		}

		if(dimension == 3) // corregir, esta mal el calculo para 3D
		{
			elongacion[0] = (1.0 + tensor_deformacion[0][0])*posicion_sin_deformar[0] + tensor_deformacion[0][1]*posicion_sin_deformar[1] + tensor_deformacion[0][2]*posicion_sin_deformar[2];
			elongacion[1] = (1.0 + tensor_deformacion[1][1])*posicion_sin_deformar[1] + tensor_deformacion[1][0]*posicion_sin_deformar[0] + tensor_deformacion[1][2]*posicion_sin_deformar[2];
			elongacion[2] = (1.0 + tensor_deformacion[2][2])*posicion_sin_deformar[2] + tensor_deformacion[2][0]*posicion_sin_deformar[0] + tensor_deformacion[2][1]*posicion_sin_deformar[1];
		}

		return elongacion;
	}
	
	
	
	
	    // Calculo de estres artificial visto en Monaghan 2000 y Gray 2001 para controlar la inestabilidad de tension
    public double[][][] Estres_Artificial(int i, int[] vecinos,   double[]  kernel, double kernel_constante, double [][][] estres_art_p_part)
    {
		int dimension  = estres_art_p_part[0].length;
        int cantidad_vecinos = vecinos.length;
        int calculo_estres_artificial = 1;
        double[][][] estres_artificial  = new double[cantidad_vecinos][dimension][dimension];

		//if(calculo_estres_artificial != 0)
		{

            int    vecino;
            double f_estres_artificial;
            double n = (double)Param.n_estres_artificial;
            double[][] matriz_temporal_A  = new double[dimension][dimension];

            for(int j=0; j< cantidad_vecinos; j++)
            {
                vecino = vecinos[j];
                matriz_temporal_A    = Util.Suma_Matrices(estres_art_p_part[i], estres_art_p_part[vecino]);
                f_estres_artificial  = Math.pow(kernel[j]/kernel_constante,n);
                estres_artificial[j] = Util.Producto_Matriz_Escalar(f_estres_artificial, matriz_temporal_A);
             }
		}

        return estres_artificial;
    }
    
    
        //calculo viscocidad artificial visto en 
    public double[][][]  Viscocidad_Artificial(int i, double[][] posicion, double[][] velocidad, double[] densidad, double[] h, int[]  vecinos, int t_actual)
    {

		
		int cantidad_vecinos = vecinos.length;
        int dimension = posicion[0].length;
		boolean calculo_visc_artificial = true;
        double[][][] viscocidad_artificial = new double[cantidad_vecinos][dimension][dimension];

		if(calculo_visc_artificial)
		{
			 double vel_sonido = Param.vel_sonido;
			int vecino;
            double valor_viscocidad = 0;
			double epsilon_cuadrado = Param.epsilon_viscosidad_artificial*Param.epsilon_viscosidad_artificial;
			double alpha  = Param.alfa_viscocidad_artificial;
			double kh_ij, c_ij, rho_ij, r, phi;	//variables auxiliares
			double beta = Param.beta_viscocidad_artificial;
			double[] v_ij = new double[dimension];
			double[] r_ij = new double[dimension];
			double k_kernel = new Kernel_y_Derivadas().K_Kernel();

			for(int j=0; j< cantidad_vecinos; j++)
			{
				vecino = vecinos[j];
				kh_ij  = (h[i]+ h[vecino])*0.5;
				c_ij   = vel_sonido;
				rho_ij = 0.5*(densidad[i] + densidad[vecino]);
				v_ij   = Util.Resta_Vectores(velocidad[i], velocidad[vecino]);
				r_ij   = Util.Resta_Vectores(posicion[i], posicion[vecino]);
				r	  = Util.Producto_Punto(r_ij,r_ij) + epsilon_cuadrado*kh_ij*kh_ij;
				phi	  = kh_ij*Util.Producto_Punto(v_ij,r_ij)/r;

				if(phi < 0)
				{
					 valor_viscocidad = (-alpha*c_ij*phi + beta*phi*phi)/(rho_ij);
				}
				if(phi >= 0)
				{
					valor_viscocidad = 0;
				}
				//viscocidad_artificial[j][0][0] = - 2.0*v_ij[0];
				//viscocidad_artificial[j][1][1] = - 2.0*v_ij[1];
                viscocidad_artificial[j] =  Util.Producto_Matriz_Escalar(valor_viscocidad, Util.Matriz_Identidad_Rango_2(dimension));
                //viscocidad_artificial[j] =  Util.Producto_Matriz_Escalar(valor_viscocidad, Util.Matriz_Componentes_Iguales(dimension, dimension, 1.0));
			}
		}
        return viscocidad_artificial;
    }
    
    	     
    public double[][] Estres_Artificial_Por_Particula(double densidad, double[][]tensor_estres)
    {
        int dimension  = tensor_estres.length;
        int calculo_estres_artificial = 1;
        double epsilon = Param.epsilon_estres_artificial;
        double inverso_densidad_al_cuadrado = 1/(densidad*densidad);
        double[][] estres_artificial  = new double[dimension][dimension];
		
        if(calculo_estres_artificial == 1)
        {
            if(tensor_estres[0][0] > 0)
            {
                estres_artificial[0][0] =  -epsilon*Math.abs(tensor_estres[0][0])*inverso_densidad_al_cuadrado;
            }
            if(tensor_estres[1][1] > 0)
            {
                estres_artificial[1][1] =  -epsilon*Math.abs(tensor_estres[1][1])*inverso_densidad_al_cuadrado;
            }
        }

        if(calculo_estres_artificial == 2)
        {
			 double angulo = 0.5*Util.Angulo_2D((tensor_estres[0][1] + tensor_estres[1][0]),(tensor_estres[0][0] - tensor_estres[1][1]));
			 double[] ejes_princ_estres_artificial = new double[dimension];
             double[] ejes_princ_estres = Estres_2D_Ejes_Principales(angulo, tensor_estres);

             if(ejes_princ_estres[0] > 0)
             {
				ejes_princ_estres_artificial[0] = - epsilon*ejes_princ_estres[0]*inverso_densidad_al_cuadrado;
             }
             if(ejes_princ_estres[1] > 0)
             {
				ejes_princ_estres_artificial[1] = - epsilon*ejes_princ_estres[1]*inverso_densidad_al_cuadrado;
             }
             
             estres_artificial = Estres_2D_Ejes_Normales(angulo, ejes_princ_estres_artificial);
        }

        return estres_artificial;
    }
    
    double Estres_Von_Mises(double estres_prin_1, double estres_prin_2)
    {
		return 	Math.sqrt(0.5*((estres_prin_1 - estres_prin_2)*(estres_prin_1 - estres_prin_2) +  estres_prin_1*estres_prin_1 + estres_prin_2*estres_prin_2   ));

	}
	
	
	double[] Tamanos_Deformacion(double[] tamanos_iniciales, double[][] tensor_deformacion)
	{
		double	coeficiente_poisson	= Param.coeficiente_poisson;
		double[] nuevos_tamanos = new double[Param.dimension + 1];
		
		nuevos_tamanos[0] = tamanos_iniciales[0]*Math.sqrt(1 + 2 *(tensor_deformacion[0][0] + tensor_deformacion[0][1]*tamanos_iniciales[1]/tamanos_iniciales[0]));
		nuevos_tamanos[1] = tamanos_iniciales[1]*Math.sqrt(1 + 2 *(tensor_deformacion[1][1] + tensor_deformacion[1][0]*tamanos_iniciales[0]/tamanos_iniciales[1]));
		nuevos_tamanos[2] = tamanos_iniciales[1]*Math.sqrt(1 + 2*coeficiente_poisson/(coeficiente_poisson - 1.0)*(tensor_deformacion[0][0] + tensor_deformacion[1][1]));
		
		return nuevos_tamanos;
	}
	
	
	double Elemento_Superficie(int borde_de_particula,  double[] vector_normal, double[] nuevos_tamanos)
	{
		double	 elemento_superficie = 0; 
		
		if(borde_de_particula == -2)  // la particula esta en el bulk
			elemento_superficie = 0;
			
		if(borde_de_particula == 0) // la particula esta en una fractura
		{
			double angulo				 = Util.Angulo_2D(vector_normal);
			double factor_angulo		 = Param.PI/180.;
			boolean limite_superior		 = angulo > 45*factor_angulo && angulo < 135*factor_angulo;
			boolean limite_inferior		 = angulo > 225*factor_angulo && angulo < 315*factor_angulo;
			
			if(limite_inferior || limite_superior) 
			/* significa que la inclinación del vector normal de la particula en la fractura tiende mas a una pared horizontal, 
			* entonces el area transversal debe ser Delta x * Delta z */
				elemento_superficie = nuevos_tamanos[0]*nuevos_tamanos[2];
			else
				elemento_superficie = nuevos_tamanos[1]*nuevos_tamanos[2];
			// Si no, entonces la pared debe estar mas inclinada a la vertical, esto es una aproximacion
		}
		
		if(borde_de_particula == 1 || borde_de_particula == 3) // paredes laterales
			elemento_superficie = nuevos_tamanos[1]*nuevos_tamanos[2];
			
		if(borde_de_particula == 2 || borde_de_particula == 4) // paredes superior e inferior
			elemento_superficie = nuevos_tamanos[0]*nuevos_tamanos[2];
		
		return	 elemento_superficie; 
	}
	
	double Energia_Superficial(double masas, double densidad, int borde_de_particula,  double[] vector_normal, double[] nuevos_tamanos)
	{
		double elemento_volumen 	= masas/densidad;
		double tension_superficial	= Param.tension_superficial;
        double elemento_superficie 	= Elemento_Superficie(borde_de_particula,  vector_normal, nuevos_tamanos)/elemento_volumen;;
        double energia_superficial  = tension_superficial*elemento_superficie;
        
        //suma la tension superficial de los laterales superior e inferior de la lamina
		energia_superficial += 2*tension_superficial*nuevos_tamanos[0]*nuevos_tamanos[1];
		
		
        
		return energia_superficial;
	}
	
	
	double Trabajo(double masas, double densidad, int borde_de_particula, double[] nuevos_tamanos, double[] traccion, double[] elongacion, double[] vector_normal)
	{
        int 	 dimension 			 = Param.dimension;   
        double	 elemento_volumen 	 = masas/densidad; 
        double	 elemento_superficie = Elemento_Superficie(borde_de_particula,  vector_normal, nuevos_tamanos)/elemento_volumen;   
		/* en realidad el elemento de superficie es elemento_superficie= elemento_volumen/separacion y el termino que se tendria es trabajo total, pero como lo que estoy comparando es 
		* trabajo por unidad de volumen, tendria que dividir entre volumen, por lo que lo hago de una vez acá */
		double[] fuerza_externa 	 = Fuerza_Externa(dimension);  
		        
		return 2*(elemento_superficie*Util.Producto_Punto(traccion, elongacion) + Util.Producto_Punto(fuerza_externa, elongacion));
	}
	
}
