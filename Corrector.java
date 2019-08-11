import java.util.List;
import java.util.ArrayList;

class Corrector
{
	
	//objetos de los metodos
	Utilidades 			Util     		= new Utilidades();
	Parametros 			Param     		= new Parametros();
	Otras_Variables 	Otras_Var 		= new Otras_Variables();
	Tasa_Variables 		Tasa_Var   		= new Tasa_Variables();
	Kernel_y_Derivadas 	Adap_y_Consis 	= new Kernel_y_Derivadas();
	Actualizar 			Actual 		  	= new Actualizar();
	Condiciones_Borde	Condicion_Borde = new Condiciones_Borde();

	
	//Variables corregidas
	double[] 	 energia_elastica ;
	double[]	 energia_cinetica;
	double[][][] tensor_deformacion;
	double[][][] tensor_estres;
	
	
	//Tasa de cambio de las variables 
	double[]     tasa_energia_elastica ;
	double[][]   velocidad;
	double[][] 	 aceleracion;

		
	//Variables relacionadas con las energias
	double 		energia_cinetica_total;
	double 		energia_elastica_total;
	double[] 	energia_auxiliar3;
	double		trabajo_total;
	double[]	trabajo;
	double		energia_superficial_total;
	double[]	energia_superficial;
	
	
	//Otras variables
	double 		delta_tiempo;
	double[][] 	estres_prin;
	double[]	nuevos_tamanos;
	
	
	//relacionadas con las fracturas y los bordes
	int[]		borde_de_particula;
	double[]	superficie;
	List<Integer> 	particulas_fractura; 
	double[][] vector_normal; 
	double[][] vector_tangente;
	List<Integer>	particulas_en_borde;
	
	
	//relacionados con la propagacion de fracturas
	double 		criterio_de_griffith;
	boolean[][] 	vecinos_hook_nuevo; 
	List<Integer>	particulas_en_borde_nuevo;
	double[]		punto_origen_fractura_nuevo;
	double[]		punto_borde_fractura_nuevo;
	double			longitud_fractura_nuevo;
	double[][]		vector_normal_nuevo;
	List<Integer>	particulas_fractura_nuevo;
	double[]		superficie_nuevo;
	int[]			borde_de_particula_nuevo;
	List<double[]>  fractura_borde_nuevo;
	List<double[]>  fractura_origen_nuevo;
	
	boolean			llego_origen_a_frontera_nuevo;
	boolean			llego_borde_a_frontera_nuevo;
		
	public Corrector
	(
		int 	 	 t_actual,
		double 	 	 delta_tiempo,
		double		 tiempo_real, 
		double 	 	 masas, 
		double[] 	 energia_elastica_anterior, 
		double[]	 tasa_energia_elastica_anterior,
		double[] 	 tasa_energia_elastica_predicha,	
		double[][] 	 posicion, 
		double[][] 	 elongacion, 
		double[][] 	 velocidad_intermedia, 
		double[] 	 h, 
		int[][] 	 vecinos,
		int[][] 	 vecinos_potenciales, 	
		int[][] 	 celda_de_particula,
		List<Integer>	particulas_en_borde,
		int[] 		 borde_de_particula, 
		double[][] 	 aceleracion_anterior,
		double[]	 separacion,
		boolean[][]  vecinos_hook,
		double[][]	 traccion,
		double[]	 densidad,
		List<Integer>	particulas_fractura,
		double[]	 punto_origen_fractura,
		double[][]	vector_normal,
		double[]	 tamano_inicial,
		double[][]	posicion_particulas_centro_de_masa,
		double[]	superficie,
		double[]	punto_borde_fractura,
		double		longitud_fractura,
		double[][]	velocidad_anterior,
		List<double[]> fractura_borde,
		List<double[]> fractura_origen,
		boolean			llego_origen_a_frontera,
		boolean			llego_borde_a_frontera
	)
	{
    	
		int N 		  = posicion.length;
		int dimension = posicion[0].length;

		boolean esta_particula_en_borde = false;
		
		// variables a integrar
		energia_elastica 	= new double[N];
		energia_cinetica 	= new double[N];
		tensor_deformacion 	= new double[N][dimension][dimension];
		tensor_estres 	  	= new double[N][dimension][dimension];
		trabajo				= new double[N];
		energia_superficial = new double[N];

		//Otras variables, internas;
		
		double[][] 	 velocidad_temporal 			= new double[N][dimension];
		double[][] 	 kernel_sin_consistencia 		= new double[N][];
		double[][][] grad_kernel_sin_consistencia 	= new double[N][1][1]; 
		double[][] 	 kernel   	 					= new double[N][];
		double[][][] grad_kernel 					= new double[N][1][1];    //se termina de inicializar en uno de los ciclos for
		double[]	 nuevas_superficies				= new double[dimension + 1];
		boolean hubo_fractura = true;
		double longitud_fractura_anterior;		
				
		//Tasa de cambio de las variables 
		velocidad     	 		= new double[N][dimension];
		aceleracion   	 		= new double[N][dimension];




		//Otras variables, externas
		energia_auxiliar3		= new double[N];
		estres_prin  			= new double[N][dimension + 1];
		criterio_de_griffith 	= new Criterio_Griffith(longitud_fractura).criterio_de_griffith;
		fractura_origen_nuevo = new ArrayList<double[]>();
		fractura_borde_nuevo = new ArrayList<double[]>();
		//criterio_de_griffith 	= Param.criterio_de_griffith_f/Math.sqrt(2.*Param.tamano[0]*Math.tan(Param.PI*longitud_fractura/(2.*Param.tamano[0])));


		for(int i = 0;  i< N; i++) // kernel y sus derivadas
		{
			Kernel_y_Derivadas kernel_y_derivadas = new Kernel_y_Derivadas(i, h , vecinos[i], posicion);
				kernel_sin_consistencia[i]      = kernel_y_derivadas.kernel;
				grad_kernel_sin_consistencia[i] = kernel_y_derivadas.gradiente;
		}

			
		for(int i = 0; i< N; i++) // kernel y sus derivadas con consistencia
		{
			kernel[i] 	   = Adap_y_Consis.Cons_Kernel(i, vecinos[i], masas, densidad, kernel_sin_consistencia[i], posicion);
			grad_kernel[i] = Adap_y_Consis.Cons_Grad(i, vecinos[i], masas, densidad, posicion, kernel_sin_consistencia[i], kernel[i], grad_kernel_sin_consistencia[i]);
		}
		
		

		//System.out.println(criterio_de_griffith);
		//criterio_de_griffith 	= 8*Math.pow(10,-3);

		
		
		for(int i = 0; i< N; i++)	//Variables necesarias para la ecuacion de movimiento
		{	
			tensor_deformacion[i] 	= Otras_Var.Deformacion_Por_Estres_Plano(i, borde_de_particula[i], vecinos[i], masas, densidad, elongacion, grad_kernel[i], vecinos_hook[i]);
			tensor_estres[i]		= Otras_Var.Tensor_Estres_A_Partir_Deformacion(i, tensor_deformacion[i]);
			estres_prin[i] 			= Util.Autovalores_Reales_Matriz_2D(tensor_estres[i]);
			estres_prin[i][2]		= Otras_Var.Estres_Von_Mises(estres_prin[i][0], estres_prin[i][1]);	
			

		}
		//System.out.println((estres_prin[1810][0] - estres_prin[1750][0]));
		
		
		//particulas_en_borde.forEach(System.out::println);
		//System.out.println( " ");
		


		
		
		vecinos_hook_nuevo 				= Util.Copiar_Matriz(vecinos_hook);
		particulas_en_borde_nuevo 		= Util.Copiar_Lista(particulas_en_borde);
		particulas_fractura_nuevo		= Util.Copiar_Lista(particulas_fractura);
		borde_de_particula_nuevo		= Util.Copiar_Vector(borde_de_particula);
		punto_origen_fractura_nuevo 	= Util.Copiar_Vector(punto_origen_fractura);
		punto_borde_fractura_nuevo		= Util.Copiar_Vector(punto_borde_fractura);
		longitud_fractura_nuevo			= longitud_fractura;
		fractura_origen_nuevo			= Util.Copiar_ListaDD(fractura_origen);
		fractura_borde_nuevo			= Util.Copiar_ListaDD(fractura_borde);
		superficie_nuevo				= Util.Copiar_Vector(superficie);
		vector_normal_nuevo				= Util.Copiar_Matriz(vector_normal);
		longitud_fractura_anterior 		= longitud_fractura;
		llego_origen_a_frontera_nuevo 	= llego_origen_a_frontera;
		llego_borde_a_frontera_nuevo 	= llego_borde_a_frontera;
		
		ArrayList<Propagacion_Fracturas> propagacion = new ArrayList<Propagacion_Fracturas>();
		
		
		

		int veces_maximas_propagacion = 0;
		
		while(hubo_fractura && tiempo_real > Param.tiempo_acustico)
		//for(int i = 0; i < 2; i++)
		{
			
			propagacion.add
			(
				new Propagacion_Fracturas
				(
					t_actual,
					vecinos, 
					estres_prin,
					posicion,
					separacion[0],
					criterio_de_griffith,
					grad_kernel,
					vecinos_hook_nuevo, 
					particulas_en_borde_nuevo,
					particulas_fractura_nuevo, 
					borde_de_particula_nuevo,
					punto_origen_fractura_nuevo,
					punto_borde_fractura_nuevo,
					longitud_fractura_nuevo,
					fractura_origen_nuevo,
					fractura_borde_nuevo,
					superficie_nuevo,
					vector_normal_nuevo,
					longitud_fractura_anterior,
					llego_origen_a_frontera_nuevo,
					llego_borde_a_frontera_nuevo 
				)
			);
			
			veces_maximas_propagacion 		= propagacion.size() - 1;
			vecinos_hook_nuevo 				= propagacion.get(veces_maximas_propagacion).vecinos_hook_nuevo;
			punto_origen_fractura_nuevo 	= propagacion.get(veces_maximas_propagacion).punto_origen_fractura_nuevo;
			punto_borde_fractura_nuevo 		= propagacion.get(veces_maximas_propagacion).punto_borde_fractura_nuevo;
			longitud_fractura_nuevo			= propagacion.get(veces_maximas_propagacion).longitud_fractura_nuevo;
			fractura_origen_nuevo			= propagacion.get(veces_maximas_propagacion).fractura_origen_nuevo;
			fractura_borde_nuevo			= propagacion.get(veces_maximas_propagacion).fractura_borde_nuevo;
			particulas_en_borde_nuevo 		= propagacion.get(veces_maximas_propagacion).particulas_en_borde_nuevo;
			vector_normal_nuevo				= propagacion.get(veces_maximas_propagacion).vector_normal_nuevo;
			borde_de_particula_nuevo		= propagacion.get(veces_maximas_propagacion).borde_de_particula_nuevo;
			superficie_nuevo				= propagacion.get(veces_maximas_propagacion).superficie_nuevo;
			hubo_fractura					= propagacion.get(veces_maximas_propagacion).hubo_fractura;
			llego_origen_a_frontera_nuevo	= propagacion.get(veces_maximas_propagacion).llego_origen_a_frontera_nuevo;
			llego_borde_a_frontera_nuevo	= propagacion.get(veces_maximas_propagacion).llego_borde_a_frontera_nuevo;
			
		}
		
		

		
		

		
		Aceleracion Ace = new Aceleracion
		(	
			particulas_en_borde_nuevo,
			borde_de_particula_nuevo, 
			t_actual,
			vecinos, 
			masas, 
			densidad,
			superficie_nuevo,
			grad_kernel,  
			tensor_estres,
			estres_prin, 
			vecinos_hook_nuevo,
			kernel,
			h,
			traccion,
			vector_normal_nuevo,
			posicion,
			separacion[0],
			velocidad_intermedia,
			velocidad_anterior
		);
			aceleracion = Ace.aceleracion;

		//Util.Imprimir_Matriz(vector_normal);
        //Util.Imprimir_Matriz(aceleracion);
   		
			
        for(int i = 0; i< N; i++)	//velocidad corregida y energias
		{
			velocidad[i] 			= Actual.Tensor_R1(velocidad_intermedia[i], aceleracion[i], 0.5*delta_tiempo);
			energia_cinetica[i]  	= 0.5*Param.densidad*Util.Producto_Punto(velocidad[i],velocidad[i]);
			energia_elastica[i]  	= 0.5*Util.Doble_Producto_Punto_Matrices(tensor_deformacion[i],tensor_estres[i] );
			nuevos_tamanos		 	= Otras_Var.Tamanos_Deformacion(tamano_inicial, tensor_deformacion[i]);
			trabajo[i]			 	= Otras_Var.Trabajo(masas, densidad[i], borde_de_particula[i], nuevos_tamanos, traccion[i],  elongacion[i], vector_normal[i]);
			energia_superficial[i] 	= Otras_Var.Energia_Superficial(masas, densidad[i], borde_de_particula[i], vector_normal[i], nuevos_tamanos);

        }   
		//System.out.println(nuevos_tamanos[1]*nuevos_tamanos[2] + " " + separacion[0]*separacion[0]);
		
		energia_elastica_total 		= Util.Suma_Componente_Vector(energia_elastica);
		energia_cinetica_total 		= Util.Suma_Componente_Vector(energia_cinetica);
		trabajo_total				= Util.Suma_Componente_Vector(trabajo);
		energia_superficial_total 	= Util.Suma_Componente_Vector(energia_superficial); 
		
		//System.out.println(trabajo_total);
		
		
		
		
		this.delta_tiempo = new Tiempo(delta_tiempo, h, velocidad, aceleracion, vecinos, t_actual).tiempo;
	}
}

