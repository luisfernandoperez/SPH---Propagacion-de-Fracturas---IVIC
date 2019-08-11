import java.util.List;

public class Tasa_Variables
{
    //este metodo solo arroja al exterior la densidad y la distancia de suavizado
	Utilidades 		Util  		= new Utilidades();
    Parametros 		Param 		= new Parametros();
	Otras_Variables Otras_Var 	= new Otras_Variables();

    
    public double Tasa_Energia_Interna2
    (
		int 		 i, 
		double	 	 masas,
		double[][] 	 velocidad, 
		double[][][] tensor_estres, 
		int[] 		 vecinos, 
		double[][] 	 grad_kernel, 
		int 		 t_actual
	)
    {
		int 		vecino 				= 0;
		int 		dimension 			= velocidad[0].length;
		double 		densidad = Param.densidad;
        double 		tasa_energia_interna 		= 0;
        double 		producto 			= 0;
		double[] 	v_ij 				= new double[dimension];
		double[][]  matriz_temporal_A   = new double[dimension][dimension];
		double[][]  producto_vel_grad  	= new double[dimension][dimension];
		
		 
		for(int j =0; j<vecinos.length; j++)
		{
			vecino 				= vecinos[j];
			matriz_temporal_A	= Util.Suma_Matrices(tensor_estres[i], tensor_estres[vecino]);
			v_ij 				= Util.Resta_Vectores(velocidad[vecino],velocidad[i]);
			producto_vel_grad 	= Util.Producto_Vectores(v_ij, grad_kernel[j]);
			producto 			= Util.Doble_Producto_Punto_Matrices(matriz_temporal_A, producto_vel_grad);
			tasa_energia_interna  	   += masas*producto/(densidad*densidad);
		}
		
		

        return tasa_energia_interna ;
    }
    

    public double Tasa_Energia_Interna
    (
		double[] 	velocidad, 
		double[] 	aceleracion, 
		double[][] 	tensor_estres, 
		double[][] 	grad_velocidades
	)
    {
		double densidad = Param.densidad;
		double[][] suma_vel = Util.Suma_Matrices(Util.Transpuesta(grad_velocidades), grad_velocidades);
       // return densidad*Util.Producto_Punto(velocidad, aceleracion) + Util.Doble_Producto_Punto_Matrices(tensor_estres, suma_vel);
        
        return - densidad*Util.Producto_Punto(velocidad, aceleracion);
    }
    public double Tasa_Energia_Helmontz
    (
		double[] 	velocidad, 
		double[] 	aceleracion, 
		double[][] 	tensor_estres, 
		double[][] 	grad_velocidades
	)
    {
		double densidad = Param.densidad;
		double[][] suma_vel = Util.Suma_Matrices(Util.Transpuesta(grad_velocidades), grad_velocidades);
        return densidad*Util.Producto_Punto(velocidad, aceleracion) + 0.5*Util.Doble_Producto_Punto_Matrices(tensor_estres, suma_vel);
        
       // return Util.Doble_Producto_Punto_Matrices(tensor_estres, grad_velocidades);
    }


}


class Aceleracion
{
	Utilidades 		Util  		= new Utilidades();
    Parametros 		Param 		= new Parametros();
	Otras_Variables Otras_Var 	= new Otras_Variables();

	double[][] aceleracion;
    
	public Aceleracion
    (	
    	List<Integer>	particulas_en_borde,
		int[] 		 	borde_de_particula, 
		int 		 	t_actual,
		int[][] 	 	vecinos, 
		double 	 	 	masas, 
		double[]	 	densidad,
 		double[] 	 	superficie,
		double[][][] 	grad_kernel,  
		double[][][] 	tensor_estres,
		double[][] 		estres_prin, 
		boolean[][]  	vecinos_hook,
		double[][]	 	kernel,
		double[]	 	h,
		double[][]	 	traccion,
		double[][]		vector_normal,
		double[][]		posicion,
		double			separacion,
		double[][]		velocidad,
		double[][]		velocidad_anterior
	)
	{
		
		boolean esta_particula_en_borde = false;

		int 		N 					= Param.N;
        int 		dimension 			= Param.dimension;
        int 		vecino				= 0;
        int 		cantidad_vecinos	= 0;
        
        double		elemento_superficie = 0; 
       	double[] 	aceleracion_hertz 	= new double[dimension];      
        double[]   	body_force     		= new double[dimension];
        double[]   	surface_force     	= new double[dimension];
        double[]   	sum_body_force		= new double[dimension];
        double[]   	sum_surface_force   = new double[dimension];
        double[]	fuerza_interna;
        double[]	traccion_i 			= new double[dimension];
        double[]	traccion_j 			= new double[dimension];
        double[]	suma_traccion		= new double[dimension];         
		double[]  	fuerza_externa 		= Otras_Var.Fuerza_Externa(dimension);
        double[][] 	matriz_temporal_A   = new double[dimension][dimension];
        double[][]  estres_i   			= new double[dimension][dimension];
		double[][]  estres_j   			= new double[dimension][dimension];

	
        aceleracion         = new double[N][dimension];
		
        for(int i = 0; i< N; i++)
		{
			cantidad_vecinos 	= vecinos[i].length;
			sum_body_force		= new double[dimension];
			sum_surface_force	= new double[dimension];

			estres_i 	= Util.Producto_Matriz_Escalar(1.0/(densidad[i]*densidad[i]), tensor_estres[i]);
			
			//if(borde_de_particula[i] == -1) System.out.println(i);
			if(superficie[i] == 0 && borde_de_particula[i]  != -2) System.out.println(i + " " + borde_de_particula[i] );
			for(int j=0; j< cantidad_vecinos; j++)
			{
				vecino 					= vecinos[i][j];
				esta_particula_en_borde = borde_de_particula[vecino] != -2;
				
				
				// 	AREA SUMATORIA VOLUMEN
				

				estres_j 				= Util.Producto_Matriz_Escalar(1.0/(densidad[vecino]*densidad[vecino]), tensor_estres[vecino]);
				matriz_temporal_A   	= Util.Suma_Matrices(estres_i, estres_j);

				// Si "vecino" es una particula perteneciente al bulk la cual interacciona via hook con "i" (incluyendo el caso cuando i=j),
				 // es decir, no hay una fractura en el medio es, entonces se hace el calculo de la fuerza con SPH

				
				
				if(vecinos_hook[i][j])
				{	 
					body_force 			= Util.Producto_Matriz_Vector(matriz_temporal_A, grad_kernel[i][j]);
					sum_body_force[0]  += masas*body_force[0] - Param.coe_friccion*(velocidad[i][0] - velocidad_anterior[i][0]);
					sum_body_force[1]  += masas*body_force[1] - Param.coe_friccion*(velocidad[i][1] - velocidad_anterior[i][1]);

				}
				//sino, no interactuan de forma cohesiva con hook, asi que se usa un potencial de esferas suaves, la ley de Hertz.
				//notese que si pasa esto, obligatoriamente las particulas estan en el borde
				else
				{
					aceleracion_hertz = Aceleracion_Hertz(esta_particula_en_borde, i, vecino, borde_de_particula[i], vecinos[i], masas, separacion, posicion);
					sum_body_force[0]  += aceleracion_hertz[0];
					sum_body_force[1]  += aceleracion_hertz[1];
				}

				

				// 	AREA SUMATORIA SUPERFICIE
				// Si la particula de prueba o su vecina estan en la superficie, se calcula la integral de superficie

				
				if(borde_de_particula[vecino] != -2)
				{
					
					//elemento_superficie	= masas*kernel[i][j]/separacion;
					elemento_superficie	= masas*kernel[i][j]/superficie[vecino];
					
					
					
					traccion_j 			= Util.Producto_Vector_Escalar(1.0/(densidad[vecino]*densidad[vecino]), traccion[vecino]);
					//traccion_i 			= Util.Producto_Vector_Escalar(1.0/(densidad[i]*densidad[i]), traccion[i]);
					traccion_i			= Componente_Traccion_i_Integral_Superficie(i, vecino, densidad[i], traccion[i], vector_normal, estres_prin[i], tensor_estres[i]);
					
					//System.out.println(i + " " + borde_de_particula[vecino] + " " + borde_de_particula[i]);
					//System.out.println(vecino + " " + vector_normal[i][1] + " " + vector_normal[vecino][1]);
					//System.out.println(traccion[i][0] + " " + traccion[i][1] + " " + densidad[vecino]);
					
					surface_force		= Util.Suma_Vectores(traccion_i, traccion_j);
					surface_force		= Util.Producto_Vector_Escalar(elemento_superficie, surface_force);
					
					
					sum_surface_force[0]   += elemento_superficie*surface_force[0];
					sum_surface_force[1]   += elemento_superficie*surface_force[1];
					
					
					//System.out.println(i + " " + vecino + " " + elemento_superficie + " " + sum_surface_force[0] + " " + sum_surface_force[1]  + " " + sum_body_force[1]  + " " + sum_body_force[1]);
				}

			}
			
		
			fuerza_interna = Util.Suma_Vectores(sum_surface_force, sum_body_force);
			aceleracion[i] = Util.Suma_Vectores(fuerza_interna, fuerza_externa);

		}
		
		//System.out.println(aceleracion[0][1] + " " +  aceleracion[4][1]);
	}
	
	
	public double[] Componente_Traccion_i_Integral_Superficie(int i, int vecino, double densidad,  double[] traccion_i_original, double[][] vector_normal, double[] estres_prin, double[][] tensor_estres)
	{
		int 		dimension 				= Param.dimension;
		double		constante_extra 		= 0;
        double		seno_rot_traccion 		= 0;
        double		seno_estres_principal 	= 0;
        double		coseno_rot_traccion 	= 0;
        double		coseno_estres_principal = 0;        
        double		angulo_estres_principal = 0;
        double		tolerancia 				= Math.pow(10, -5);  
		double[]	traccion_i;
		double[]	vector_extra;
		double[]	vector_normal_rotado;
		double[][]	rotacion_traccion;
		double[][]  rotacion_estres; 
		double[][]	auxiliar_rotacion = new double[][] { {0, -1}, {1, 0} };
		
		
		seno_rot_traccion 	= - Util.Producto_Cruz(vector_normal[vecino], vector_normal[i])[0];
		coseno_rot_traccion = Util.Producto_Punto(vector_normal[i], vector_normal[vecino]);
		rotacion_traccion 	= Matriz_Rotacion(seno_rot_traccion, coseno_rot_traccion);
		
		
		traccion_i 			= Util.Producto_Matriz_Vector(rotacion_traccion, traccion_i_original);
		
		if(estres_prin[0] != estres_prin[1])
			angulo_estres_principal = 0.5*Math.asin((tensor_estres[0][1] + tensor_estres[1][0])/(estres_prin[0] - estres_prin[1]));
		else
			angulo_estres_principal = 0;
			
		seno_estres_principal	= Math.sin(angulo_estres_principal);
		coseno_estres_principal = Math.cos(angulo_estres_principal);
		rotacion_estres			= Matriz_Rotacion(seno_estres_principal, coseno_estres_principal);
		vector_normal_rotado	= Util.Producto_Matriz_Vector(rotacion_estres, vector_normal[i]);
		
		if(vector_normal_rotado[0] > tolerancia )
			constante_extra  = traccion_i_original[0]/vector_normal_rotado[0]; 
		
		if(vector_normal_rotado[1] > tolerancia )	
			constante_extra += traccion_i_original[1]/vector_normal_rotado[1];
			 
		constante_extra = seno_rot_traccion*constante_extra;
		
		vector_extra 	= Util.Producto_Matriz_Vector(auxiliar_rotacion, vector_normal[i]);
		vector_extra 	= Util.Producto_Vector_Escalar(constante_extra, vector_extra);
		
		traccion_i		= Util.Resta_Vectores(traccion_i, vector_extra);
		
		return Util.Producto_Vector_Escalar(1.0/(densidad*densidad),traccion_i);
	}
	
	public double[][] Matriz_Rotacion(double seno, double coseno)
	{
		int 		dimension 	 = Param.dimension;
		double[][]	matriz_rotacion = new double[dimension][dimension];

		matriz_rotacion[0][0]	= coseno;
		matriz_rotacion[0][1]	= -seno;
		matriz_rotacion[1][0]	= seno;
		matriz_rotacion[1][1]	= coseno;
		
		return matriz_rotacion;
	}
	
	public double[] Aceleracion_Hertz
    (
    	boolean 	esta_particula_en_borde,
		int 		i,
		int			j, 
		int 		borde_de_particula, 
		int[] 		vecinos,
		double 		masas,
		double 	separacion,
		double[][] 	posicion
    )
    {
		int 		dimension 			= Param.dimension;
        int 		vecino;
        int 		cantidad_vecinos 	= vecinos.length;
        double		dx;
        double		dy;
        double		angulo				= 0;
        double		modulo_aceleracion	= 0;
        double 		distancia;
        double		deformacion			= 0;		//esta variable NO es la misma de la variable elongacion de los metodos anteriores
        double 		densidad 			= Param.densidad;
		double 		radio 				= 0.5*separacion;
		double		radio_efectivo 		= 0.5*radio;	//el radio efectivo es 1/reff = 1/r1 + 1/r2, pero como la separacion entre las particulas es igual...
		double	 modulo_compresibilidad = Param.modulo_compresibilidad;
		double 		modulo_rigidez		= Param.modulo_rigidez;
		double		cons_elasticidad	= Param.cons_hertz;
		double		cons_kappa			= 2*Math.sqrt(radio_efectivo)/3*cons_elasticidad/masas;
		double		tolerancia_defom	= Math.pow(10,-15);
        double[]   	aceleracion         = new double[dimension];



			if(i != j)
			{
				dx 					= posicion[i][0] - posicion[j][0];
				dy 					= posicion[i][1] - posicion[j][1];
				distancia 			= Math.sqrt(dx*dx + dy*dy);
				deformacion 		= 2*radio - distancia;
				if(Math.abs(deformacion) < tolerancia_defom)
					deformacion = 0;
				if(deformacion > 0)		
					modulo_aceleracion 	=   cons_kappa*Math.pow(deformacion, 1.5);
				else
					modulo_aceleracion 	= 0;
				angulo 				= Util.Angulo_2D(dx, dy);
				aceleracion[0] 	   += modulo_aceleracion*Math.cos(angulo);
				aceleracion[1] 	   += modulo_aceleracion*Math.sin(angulo);
				//if(i == 0 && j == 1) System.out.println(j+ " " + aceleracion[0]  + " " + dx + " " + deformacion );
			}			
			
        
		return  aceleracion;
	}


}	
