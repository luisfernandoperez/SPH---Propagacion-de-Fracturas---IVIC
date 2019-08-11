public class Condiciones_Borde
{

    Utilidades 		Util 		= new Utilidades();
    Parametros 		Param 		= new Parametros();
    Otras_Variables Otras_Var 	= new Otras_Variables();
	Actualizar 		Actual 		= new Actualizar();


	
	public double[] Condiciones_Desplazamiento(int i, boolean esta_particula_en_borde, int borde_de_particula, double[] elongacion)
    {
		//para condiciones_borde: borde izquierdo = 0, borde superior = 1, borde derecho = 2, borde inferior = 3
		int[] 	 	tipo_cond_borde 		= Param.tipo_cond_borde;
		double[] 	desplazamiento  		= Util.Copiar_Vector(elongacion);
		double[][] 	vector_desplazamiento 	= Param.vector_desplazamiento;

		if(esta_particula_en_borde && borde_de_particula > -1)
		{
			if(tipo_cond_borde[borde_de_particula] == 3 || tipo_cond_borde[borde_de_particula] == 4 )
			{
				if(borde_de_particula%2 == 0)
				{
					desplazamiento[0] = vector_desplazamiento[borde_de_particula][0];
					desplazamiento[1] = vector_desplazamiento[borde_de_particula][1];

				}
		
				if((borde_de_particula+1)%2 == 0)
				{
					desplazamiento[0] = vector_desplazamiento[borde_de_particula][0];
					desplazamiento[1] = vector_desplazamiento[borde_de_particula][1];
				}
			}
		}
		return desplazamiento;
	}

  
  
	  

	 public double[][] Condiciones_Traccion(int tiempo, int i, boolean esta_particula_en_borde, int borde_de_particula, double[] vector_normal,  double[] posicion_particulas_centro_de_masa,  double[][] tensor_deformacion)
     {	
		boolean 	 elastoestatica 			= Param.elastoestatica;
		int[] 		 tipo_cond_borde 			= Param.tipo_cond_borde;
		int[] 		 orden_condiciones_borde 	= Param.orden_condiciones_borde;
		double[][]   vector_traccion			= Establecer_Vector_Traccion(posicion_particulas_centro_de_masa);
		double[][] 	 tensor_estres 	  			= Otras_Var.Tensor_Estres_A_Partir_Deformacion(i, tensor_deformacion);

		/* La clasificacion de los bordes es la siguiente:
		 * -2 para particulas del bulk - que no estan en el borde-
		 *  0 para particulas que esten en una fractura
		 *  1 para pared izquierda
		 *  2 para pared superior
		 *  3 para pared derecha
		 *  4 para pared inferior
		 *  -1 para particulas a las que no se les haya podido encontrar pared
		 */
		
		// solo se le aplican condiciones de borde a aquellas que tengan un borde definido 
		if(borde_de_particula > -1) 
		{
			// solo se le aplican condiciones de borde de traccion a aquellas que en Parametros.java se les haya puesto
			if(tipo_cond_borde[borde_de_particula] == 1 || tipo_cond_borde[borde_de_particula] == 2)
			{
				tensor_estres = Estres_Borde(i, borde_de_particula, vector_normal,  vector_traccion, tensor_estres);
			}
		}
		//if(borde_de_particula == -1) 
		//		tensor_estres = new double[dimension][dimension];
		
		
			
		return tensor_estres;
	}

	 public double[][] Condiciones_Traccion(int i, int borde_de_particula, double tiempo, double[] vector_normal,  double[] posicion_particulas_centro_de_masa, double[][] tensor_estres_anterior, double[][] tasa_estres_anterior, double[][] tasa_estres_predicha)
     {	
		int[] 		 tipo_cond_borde 			= Param.tipo_cond_borde;
		int[] 		 orden_condiciones_borde 	= Param.orden_condiciones_borde;
		double[][]   vector_traccion			= Establecer_Vector_Traccion(posicion_particulas_centro_de_masa);
		double[][] 	 tensor_estres 	  			= Actual.Tensor_R2_Final_PC(tensor_estres_anterior, tasa_estres_predicha, tasa_estres_anterior, tiempo);


		// solo se le aplican condiciones de borde a aquellas que tengan un borde definido 
		if(borde_de_particula > -1) 
		{
			// solo se le aplican condiciones de borde de traccion a aquellas que en Parametros.java se les haya puesto
			if(tipo_cond_borde[borde_de_particula] == 1 || tipo_cond_borde[borde_de_particula] == 2)	
				tensor_estres = Estres_Borde(i, borde_de_particula, vector_normal,  vector_traccion, tensor_estres);
		}
		
			
		return tensor_estres;
	}
	
	public double[][] Estres_Borde(int i, int borde_de_particula, double[] vector_normal, double[][] vector_traccion, double[][] tensor_estres )
    {
		//para condiciones_borde: borde izquierdo = 0, borde superior = 1, borde derecho = 2, borde inferior = 3
		double tolerancia = Math.pow(10,-5);

		boolean nx_cero = Math.abs(vector_normal[0]) <= tolerancia;
		boolean ny_cero = Math.abs(vector_normal[1]) <= tolerancia;
		int dimension = vector_traccion[0].length;
		double[][] estres_borde = new double[dimension][dimension];
		
		int[] 	tipo_cond_borde = Param.tipo_cond_borde;
	
		/*
		if((borde_de_particula + 1)%2 == 0)	// paredes impares son laterales
		{
			estres_borde[0][0] =  tensor_estres[0][0] +  vector_traccion[borde_de_particula][0];
			estres_borde[0][1] =  tensor_estres[0][1] +   vector_traccion[borde_de_particula][1];
			estres_borde[1][0] =   estres_borde[0][1];
			estres_borde[1][1]	= tensor_estres[1][1] ;
		}

		if(borde_de_particula%2 == 0)	//paredes pares son superior o inferior
		{
			estres_borde[1][0] = tensor_estres[1][0] +  vector_traccion[borde_de_particula][0];
			estres_borde[0][1] = estres_borde[1][0];
			estres_borde[1][1] = tensor_estres[1][1] +  vector_traccion[borde_de_particula][1];
			estres_borde[0][0] = tensor_estres[0][0];
		}
		
		
		
		if(borde_de_particula == 0)	// paredes de fractura, que en esta construccion son inferior/superior
		{
			estres_borde[1][0] = 0;
			estres_borde[0][1] = 0;
			estres_borde[1][1] = 0;
			estres_borde[0][0] = 0;
		}
		*/
	
		
		if((borde_de_particula + 1)%2 == 0)	// paredes impares son laterales
		{
			estres_borde[0][0] = vector_traccion[borde_de_particula][0];
			estres_borde[0][1] = vector_traccion[borde_de_particula][1];
			estres_borde[1][0] =   estres_borde[0][1];
			estres_borde[1][1]	= 0 ;
		}

		if(borde_de_particula != 0 && borde_de_particula%2 == 0)	//paredes pares son superior o inferior
		{
			estres_borde[1][0] = vector_traccion[borde_de_particula][0];
			estres_borde[0][1] = estres_borde[1][0];
			estres_borde[1][1] = vector_traccion[borde_de_particula][1];
			estres_borde[0][0] = 0;
		}
		
		
		
		if(borde_de_particula == 0)	// paredes de fractura, que en esta construccion son inferior/superior
		{
			//estres_borde[1][0] = vector_traccion[borde_de_particula][0];
			//estres_borde[0][1] = estres_borde[1][0];
			//estres_borde[1][1] = vector_traccion[borde_de_particula][1];
			//estres_borde[0][0] = tensor_estres[0][0];
			estres_borde[1][0] = 0;
			estres_borde[0][1] = 0;
			estres_borde[1][1] = 0;
			estres_borde[0][0] = 0;
		}
		
		
		
		return estres_borde;	
	}
	


	public double[][] Vector_Traccion(int[] borde_de_particula, double[][] posicion_particulas_centro_de_masa)
	{
		int N 						= Param.N;
		int dimension 				= Param.dimension;
		double[][] vector_traccion ;
		double[][] traccion 		= new double[N][dimension];
		
		/* La clasificacion de los bordes es la siguiente:
		 * -2 para particulas del bulk - que no estan en el borde-
		 *  0 para particulas que esten en una fractura
		 *  1 para pared izquierda
		 *  2 para pared superior
		 *  3 para pared derecha
		 *  4 para pared inferior
		 *  -1 para particulas a las que no se les haya podido encontrar pared
		 */
		 
		for(int i =0; i < N; i++)
		{
			vector_traccion = Establecer_Vector_Traccion(posicion_particulas_centro_de_masa[i]);
			
			if(borde_de_particula[i] > -1)
				traccion[i] = Util.Copiar_Vector(vector_traccion[borde_de_particula[i]]);
				
		}
		
		//Util.Imprimir_Matriz(traccion);
		
		return traccion; 
	}
	

	public double[][] Establecer_Vector_Traccion(double[] posicion_particulas_centro_de_masa)
    {

		//SE SUPONE QUE EN UNA DIMENSION LA BARRA ESTA ACOSTADA EN EL EJE X
		//la primera componente indica sobre cual pared esta actuando, borde izquierdo = 0, borde superior = 1, borde derecho = 2, borde inferior = 3
		//la segunda indica las componentes espaciales, 0 para "x" y 1 para "y"
		int dimension 					= Param.dimension;
		int numero_condiciones_borde 	= Param.numero_condiciones_borde;
		double[][] vector_traccion 		= new double[numero_condiciones_borde][dimension];
		
		double	tension_critica_dimensional = Param.tension_critica_dimensional;
		
		//System.out.println(tension_critica_dimensional + " " + tension_critica);
		
		for(int condiciones_borde = 0; condiciones_borde < numero_condiciones_borde; condiciones_borde++)
		{
			// borde de fracturas = 0
			if(condiciones_borde == 0)
			{
				vector_traccion[condiciones_borde][0] = 0*Math.pow(10,5);
				vector_traccion[condiciones_borde][1] = 0*Math.pow(10,5);
			}

			// borde izquierdo  = 1
			if(condiciones_borde == 1)
			{
				vector_traccion[condiciones_borde][0] = -0*tension_critica_dimensional;
				vector_traccion[condiciones_borde][1] = 0*Math.pow(10,3);
			}
			
			// borde superior = 2		
			if(condiciones_borde == 2)
			{
				vector_traccion[condiciones_borde][0] = 0*Math.pow(10,5);
				vector_traccion[condiciones_borde][1] = -tension_critica_dimensional;
			}
			
			//borde derecho = 3
			if(condiciones_borde == 3)
			{
				vector_traccion[condiciones_borde][0] = 0*tension_critica_dimensional;
				vector_traccion[condiciones_borde][1] = 0*Math.pow(10,3);
			}
			
			//borde inferior = 4
			if(condiciones_borde == 4)
			{
				vector_traccion[condiciones_borde][0] = 0*Math.pow(10,5);
				vector_traccion[condiciones_borde][1] = tension_critica_dimensional;
			}

			//condiciones de borde para fracturas es
			vector_traccion[condiciones_borde] = Util.Producto_Vector_Escalar(1/Param.modulo_compresibi_real,vector_traccion[condiciones_borde]);
			
		}

		
		return vector_traccion;	
	}
	
	
}
