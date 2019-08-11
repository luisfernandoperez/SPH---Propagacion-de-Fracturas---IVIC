import java.io.*;
import java.util.*;


public class Imprimir
{
	Utilidades Util  = new Utilidades();
	Parametros Param = new Parametros();


	boolean condicion_imprimir = Param.condicion_imprimir;
	// Variables de impresion
	String[]   	nombres_serie_temporal;
	String[]   	nombres_serie_espacial;
	int 		N 		  						= Param.N;
	int 		particula_prueba 				= Param.particula_prueba;
	int        	numero_variables_guardar 	   	= Param.numero_variables_guardar;
	int        	variables_extra_serie_temporal  = Param.variables_extra_serie_temporal;
	int 		imprimir_cada_x_pasos			= Param.imprimir_cada_x_pasos;
	double[]   	serie_temporal;
	double[][] 	serie_espacial;
	
	public Imprimir
	(
		int 		 t_actual,
		int 		 tmax,
		double 		 tiempo_real,
		double[][] 	 posicion,
		double[][] 	 elongacion,	  	  
		double[] 	 posicion_centro_de_masa,  	  
		double[][] 	 posicion_particulas_centro_de_masa,
		double[][] 	 velocidad,
		double[] 	 velocidad_centro_de_masa,	  
		double[][] 	 velocidad_particulas_centro_de_masa,	    
		double[][] 	 aceleracion,  
		double[] 	 energia_elastica ,
		double[] 	 tasa_energia_elastica ,		   
		double[][][] tensor_deformacion,	  
		double[][][] tensor_estres,  	  
		double[][] 	 estres_prin,	  		
		double[]  	 h,				
		double		energia_cinetica_total,
		double		energia_elastica_total,
		double[] 	energia_auxiliar3,
		boolean[][] vecino_hook,
		double[]	densidad,
		double[] 	punto_origen_fractura,
		double		separacion,
		double		trabajo_total,
		double[]	trabajo,
		double		energia_superficial_total,
		double[]	energia_superficial,
		double[]	energia_cinetica,
		double		criterio_de_griffith,
		double		longitud_fractura,
		double		longitud_fractura_anterior
	)
	{
		 
		if(t_actual <= tmax || tmax < 0)	// Area de "imprersion" de las magnitudes calculadas en el leapfrog
		{
			
			//System.out.println(particula_prueba + " " + N);
			
			// IMPRESION EN ARCHIVOS .TXT
			
			//Serie Temporal
			nombres_serie_temporal = Nombre_Serie_Temporal();
			serie_temporal 	       = Serie_Temporal
									(
										particula_prueba, 
										tiempo_real, 
										posicion, 
										elongacion,  
										posicion_centro_de_masa, 
										posicion_particulas_centro_de_masa, 
										velocidad,  
										velocidad_centro_de_masa, 
										velocidad_particulas_centro_de_masa,  
										aceleracion,  
										energia_elastica , 
										tensor_deformacion, 
										tensor_estres, 									  
										estres_prin, 
										h, 
										energia_cinetica_total, 
										energia_elastica_total, 
										energia_auxiliar3,
										densidad,
										separacion,
										trabajo_total,
										trabajo,
										energia_superficial_total,
										energia_superficial,
										energia_cinetica
									);

			
			if(t_actual == 0)
			{
				String directorio_serie_temporal = Param.directorio_serie_temporal;
				
				//File f = new File("Main_SPH.java");
				//System.out.println(f.getParentFile().getName());
					//	System.out.println( new String(Param.directorio_superior + "serie_temporal/") + " " + Param.directorio_superior);

				
				new File(directorio_serie_temporal).mkdirs();	//Crea la carpeta
				
				try 	//crea el archivo y borra los viejos si ya existen otros
				{ 
					String serie_temporal 		     = Param.serie_temporal;	
					FileOutputStream 	crear_tp 	 = new FileOutputStream(serie_temporal);
					OutputStreamWriter 	salida_tp 	 = new OutputStreamWriter(crear_tp);    
					Writer 				escribir_tp  = new BufferedWriter(salida_tp);
					PrintWriter 		cabecera_tp  = new PrintWriter(salida_tp);
					escribir_tp.write(" ");

					for(int entradas = 0; entradas < nombres_serie_temporal.length; entradas++)
					{
						cabecera_tp.write(" " + nombres_serie_temporal[entradas]);
					}
					cabecera_tp.write(System.getProperty("line.separator"));
					escribir_tp.close();
				} 
				catch (IOException e) { System.err.println("No se han creado los archivos de la serie temporal");}
							
	
			}
			
			Util.Escribir_en_Tiempo(serie_temporal, nombres_serie_temporal);
			
			
		
			//Serie espacial, se imprime cada multiplo de imprimir_cada_x_pasos			
			if(t_actual%imprimir_cada_x_pasos == 0)
			{
				nombres_serie_espacial = Nombre_Serie_Espacial();
				serie_espacial	       = Serie_Espacial
										(
											posicion, 
											elongacion, 
											posicion_particulas_centro_de_masa, 
											velocidad,  
											velocidad_particulas_centro_de_masa,  
											aceleracion,    
											energia_elastica , 
											tasa_energia_elastica ,  
											tensor_deformacion, 
											tensor_estres,    
											estres_prin,     
											h
										);

				String directorio_serie_espacial = Param.directorio_serie_espacial;
				
				File carpeta_serie_espacial = new File(directorio_serie_espacial);
							
				carpeta_serie_espacial.mkdirs();	//Crea la carpeta	
			
				try 	//crea el archivo y borra los viejos si ya existen otros
				{ 
					if(t_actual == 0)
						Util.Limpiar_Archivos(carpeta_serie_espacial, "txt");
					
					String serie_espacial 		     = Param.serie_espacial + "_" + String.valueOf(t_actual) + ".txt";

					FileOutputStream 	crear_ep 	 = new FileOutputStream(serie_espacial);
					OutputStreamWriter 	salida_ep 	 = new OutputStreamWriter(crear_ep);    
					Writer 				escribir_ep  = new BufferedWriter(salida_ep);
					PrintWriter 		cabecera_ep  = new PrintWriter(salida_ep);
					escribir_ep.write(" ");

					for(int entradas = 0; entradas < nombres_serie_espacial.length; entradas++)
					{
						cabecera_ep.write(" " + nombres_serie_espacial[entradas]);
					}
					cabecera_ep.write(System.getProperty("line.separator"));
					escribir_ep.close();
				} 
				catch (IOException e) 
				
				{ System.err.println("No se han creado los archivos de la serie espacial");}

				//System.out.println(serie_espacial.length);
				Util.Escribir_en_Espacio(serie_espacial, nombres_serie_espacial, t_actual);
			}		 
			
			
			
			
			// IMPRESION EN LA TERMINAL					

			if(condicion_imprimir)
			{
				System.out.println
				(
				t_actual  + " " +
				tiempo_real/Param.tiempo_acustico 
				//h[0]
				//Util.Error_Medio_Cuadratico(funcion_SPH, funcion) + " " +
				//numero_vecinos[0] + " " +	
				//tensor_estres[0][0][0]
				//tensor_deformacion[0][0][0]
				//energia_elastica _helmontz[0]		
				//suma_kernel[30]
				//for(int i = 0; i< N; i++){
				//velocidad[0][0] + " " + velocidad[N-1][0]
				//energia_elastica [0] + " " + energia_elastica [(Param.N-(int)Math.sqrt(Param.N))/2]
				//posicion[particula_prueba][1]
				//posicion[923][1] + " " + posicion[1019][1] 
				//vecinos[0].length + " " + vecinos[particula_medio[0]].length
                // kernel[927][0] + " " + vecinos[927][0]
				//densidad[0] + " " + densidad[N-1]
				//energia_elastica [0] + " " + energia_elastica [medio] + " " + energia_elastica [N-1] 
				// posicion[0][0] + " " +  posicion[0][1]
				//aceleracion[0][0] + " " + aceleracion[N-1][0] 
				//elongacion[0][0] + " " + elongacion[1][0] 
				//vel_sonido[923] + " " +  vel_sonido[1019]
				//presion[923] + " " +  presion[1019]
				//rmse = Util.Error(rmse, promedio,promedio);
				//rmse + " " + promedio
				 + " " + tiempo_real
				);
			}
			
			
			if(t_actual == 0)
			{
				if(Param.obstaculo == -1)
				{
					System.out.println(" ");
					System.out.println("Las coordenadas del punto de fractura son:");
					Util.Imprimir_Vector(punto_origen_fractura);
				}
				
				if(Param.Imprimir_Estres_Maximo_Y_Griffith)
				{
					System.out.println(" ");
					System.out.println("La Tracción externa que indice sobre las paredes del medio es: " + Param.tension_critica_dimensional + " Pa o " + Param.tension_critica*Param.factor_estres + " adimensional");
					System.out.println("La condición inicial de Griffith es: " + criterio_de_griffith*Param.modulo_compresibi_real + " Pa o " + criterio_de_griffith + " adimensional");
				}
			}
			
			if(longitud_fractura_anterior != longitud_fractura && Param.Imprimir_Info_Fractura)
					System.out.println(longitud_fractura + " " +  tiempo_real + " " + t_actual + " " + criterio_de_griffith);
					
					
		}
	}
	
	public double[] Serie_Temporal
	(
		int 		 particula, 		
		double 		 tiempo_real,
		double[][] 	 posicion,
		double[][] 	 elongacion,	  	  
		double[] 	 posicion_centro_de_masa,  	  
		double[][] 	 posicion_particulas_centro_de_masa,
		double[][] 	 velocidad,
		double[] 	 velocidad_centro_de_masa,	  
		double[][] 	 velocidad_particulas_centro_de_masa,	   
		double[][] 	 aceleracion,  
		double[] 	 energia_elastica ,  
		double[][][] tensor_deformacion,	  
		double[][][] tensor_estres,  
		double[][] 	 estres_prin,	  		
		double[]  	 h,
		double		energia_cinetica_total,
		double		energia_elastica_total,
		double[] 	energia_auxiliar3,
		double[] 	densidad,
		double		separacion,
		double		trabajo_total,
		double[]	trabajo,
		double		energia_superficial_total,
		double[]	energia_superficial,
		double[]	energia_cinetica
	)
	{
		double[] serie_temporal	= new double[numero_variables_guardar + variables_extra_serie_temporal];
		int par_medio_iz = Param.par_medio_iz;
		int par_medio_de = Param.par_medio_de;

		double[][] trans_pos = Util.Transpuesta(posicion);

		serie_temporal[0]  = tiempo_real;
		serie_temporal[1]  = posicion[particula][0];
		serie_temporal[2]  = posicion[particula][1];
		serie_temporal[3]  = elongacion[particula][0];		  
		serie_temporal[4]  = elongacion[particula][1];
		serie_temporal[5]  = Util.Maximo(trans_pos[0]) - Util.Minimo(trans_pos[0]);
		serie_temporal[6]  = Util.Maximo(trans_pos[1]) - Util.Minimo(trans_pos[1]);  
		serie_temporal[7]  = posicion_centro_de_masa[0];		  
		serie_temporal[8]  = posicion_centro_de_masa[1];		  
		serie_temporal[9]  = posicion_particulas_centro_de_masa[particula][0];
		serie_temporal[10] = posicion_particulas_centro_de_masa[particula][1];		  
		serie_temporal[11] = velocidad[particula][0];
		serie_temporal[12] = velocidad[particula][1];
		serie_temporal[13] = velocidad_centro_de_masa[0];		  
		serie_temporal[14] = velocidad_centro_de_masa[1];
		serie_temporal[15] = velocidad_particulas_centro_de_masa[particula][0];
		serie_temporal[16] = velocidad_particulas_centro_de_masa[particula][1];		  
		serie_temporal[17] = 0;		  
		serie_temporal[18] = aceleracion[particula][0];		  
		serie_temporal[19] = aceleracion[particula][1];
		serie_temporal[20] = densidad[particula];	
		serie_temporal[21] = Param.tamano[0];
		serie_temporal[22] = Param.tamano[1];
		serie_temporal[23] = separacion;		  
		serie_temporal[24] = tensor_deformacion[particula][0][0];
		serie_temporal[25] = tensor_deformacion[particula][0][1];		  
		serie_temporal[26] = tensor_deformacion[particula][1][1];  		  
		serie_temporal[27] = tensor_estres[particula][0][0];
		serie_temporal[28] = tensor_estres[particula][0][1];
		serie_temporal[29] = tensor_estres[particula][1][1];
		serie_temporal[30] = energia_cinetica[particula];	  
		serie_temporal[31] = energia_elastica[particula];	  
		serie_temporal[32] = trabajo[particula];
		serie_temporal[33] = energia_superficial[particula];
		serie_temporal[34] = energia_cinetica_total;
		serie_temporal[35] = energia_elastica_total;
		serie_temporal[36] = trabajo_total;
		serie_temporal[37] = energia_superficial_total; 
		serie_temporal[38] = 0;
		serie_temporal[39] = energia_auxiliar3[particula];	  
		serie_temporal[40] = estres_prin[particula][0];		  
		serie_temporal[41] = estres_prin[particula][1];		  
		serie_temporal[42] = estres_prin[particula][2];
		serie_temporal[43] = 0;	
		serie_temporal[44] = 0;
		serie_temporal[45] = 0;
		serie_temporal[46] = 0;		  
		serie_temporal[47] = 0;
		serie_temporal[48] = 0;
		serie_temporal[49] = -1;
		serie_temporal[50] = h[particula];		

		return serie_temporal;
	}

	public String[] Nombre_Serie_Temporal()
	{
		String[] nombres = new String[numero_variables_guardar + variables_extra_serie_temporal];
		
		nombres[0]  = "tiempo";		
		nombres[1]  = "posicion.en.x";			
		nombres[2]  = "posicion.en.y";
		nombres[3]  = "elongacion.en.x";
		nombres[4]  = "elongacion.en.y";            
		nombres[5]  = "tamanio.en.x";            
		nombres[6]  = "tamanio.en.y";          
		nombres[7]  = "posicion.centro.de.masa.x";
		nombres[8]  = "posicion.centro.de.masa.y";         
		nombres[9]  = "posicion.particulas.con.respecto.centro.de.masa.x";         
		nombres[10] = "posicion.particulas.con.respecto.centro.de.masa.y";   	
		nombres[11] = "velocidad.en.x";			
		nombres[12] = "velocidad.en.y";
		nombres[13] = "velocidad.centro.de.masa.x";
		nombres[14] = "velocidad.centro.de.masa.y";          
		nombres[15] = "velocidad.particulas.con.respecto.centro.de.masa.x";       
		nombres[16] = "velocidad.particulas.con.respecto.centro.de.masa.y";         
		nombres[17] = "angulo.centro.de.masa";
		nombres[18] = "aceleracion.en.x";            
		nombres[19] = "aceleracion.en.y"; 
		nombres[20] = "densidad";		
		nombres[21] = "tasa.densidad";			
		nombres[22] = "energia_elastica .helmontz";
		nombres[23] = "tasa.energia_elastica .helmontz";		
		nombres[24] = "energia_elastica .interna";
		nombres[25] = "tasa.energia_elastica .interna";	
		nombres[26] = "energia_elastica .cinetica";
		nombres[27] = "deformacion.x.x";            
		nombres[28] = "deformacion.x.y";            
		nombres[29] = "deformacion.y.y";          
		nombres[30] = "estres.x.x";
		nombres[31] = "estres.x.y";            
		nombres[32] = "estres.y.y"; 
		nombres[33] = "energia_cinetica_total";
		nombres[34] = "energia_elastica_total";		
		nombres[35] = "energia_elastica_total";
		nombres[36] = "energia_cinetica_total"; 	
		nombres[37] = "energia_auxiliar3";	            
		nombres[38] = "sin.uso";        
		nombres[39] = "angulo.estres.principal";   
		nombres[40] = "estres.principal.primero"; 
		nombres[41] = "estres.principal.segundo";                
		nombres[42] = "estres.von.mises";     
		nombres[43] = "suma.kernel";		
		nombres[44] = "suma.kernel.con.consistencia";			
		nombres[45] = "suma.gradiente.kernel.x";
		nombres[46] = "suma.gradiente.kernel.y";
		nombres[47] = "suma.gradiente.kernel.x.con.consistencia";           
		nombres[48] = "suma.gradiente.kernel.y.con.consistencia";    
		nombres[49] = "sin.uso";           
		nombres[50] = "distancia.suavizado";   

		return nombres;
	}

	public double[][] Serie_Espacial
	(
		double[][] 	 posicion,
		double[][] 	 elongacion,	  	    
		double[][] 	 posicion_particulas_centro_de_masa,
		double[][] 	 velocidad,  
		double[][] 	 velocidad_particulas_centro_de_masa,	  	  
		double[][] 	 aceleracion,  
		double[] 	 energia_elastica ,
		double[] 	 tasa_energia_elastica ,		   
		double[][][] tensor_deformacion,	  
		double[][][] tensor_estres,  
		double[][] 	 estres_prin,
		double[]  	 h
	)
	{
		double[][] serie_espacial	= new double[numero_variables_guardar][N];
				
		serie_espacial[0]  = Util.Copiar_Vector(Util.Transpuesta(posicion)[0]);
		serie_espacial[1]  = Util.Copiar_Vector(Util.Transpuesta(posicion)[1]);
		serie_espacial[2]  = Util.Copiar_Vector(Util.Transpuesta(elongacion)[0]);		  
		serie_espacial[3]  = Util.Copiar_Vector(Util.Transpuesta(elongacion)[1]);	    
		serie_espacial[4]  = Util.Copiar_Vector(Util.Transpuesta(posicion_particulas_centro_de_masa)[0]);
		serie_espacial[5]  = Util.Copiar_Vector(Util.Transpuesta(posicion_particulas_centro_de_masa)[1]);		  
		serie_espacial[6]  = Util.Copiar_Vector(Util.Transpuesta(velocidad)[0]);
		serie_espacial[7]  = Util.Copiar_Vector(Util.Transpuesta(velocidad)[1]);
		serie_espacial[8]  = Util.Copiar_Vector(Util.Transpuesta(velocidad_particulas_centro_de_masa)[0]);
		serie_espacial[9]  = Util.Copiar_Vector(Util.Transpuesta(velocidad_particulas_centro_de_masa)[1]);		    
		serie_espacial[10] = Util.Copiar_Vector(Util.Transpuesta(aceleracion)[0]);		  
		serie_espacial[11] = Util.Copiar_Vector(Util.Transpuesta(aceleracion)[1]);
		serie_espacial[12] = Util.Vector_Componentes_Iguales(N, Param.densidad);	
		serie_espacial[13] = Util.Vector_Componentes_Iguales(N,0.0);
		serie_espacial[14] = Util.Copiar_Vector(energia_elastica );
		serie_espacial[15] = Util.Copiar_Vector(tasa_energia_elastica );		  
		serie_espacial[16] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[17] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[18] = Util.Vector_Componentes_Iguales(N,0.0);	 
		serie_espacial[19] = Util.Copiar_Vector(Util.Transpuesta_Primer_Indice_Al_Final(tensor_deformacion)[0][0]);
		serie_espacial[20] = Util.Copiar_Vector(Util.Transpuesta_Primer_Indice_Al_Final(tensor_deformacion)[0][1]);
		serie_espacial[21] = Util.Copiar_Vector(Util.Transpuesta_Primer_Indice_Al_Final(tensor_deformacion)[1][1]);  
		serie_espacial[22] = Util.Copiar_Vector(Util.Transpuesta_Primer_Indice_Al_Final(tensor_estres)[0][0]);	  
		serie_espacial[23] = Util.Copiar_Vector(Util.Transpuesta_Primer_Indice_Al_Final(tensor_estres)[0][1]);	  
		serie_espacial[24] = Util.Copiar_Vector(Util.Transpuesta_Primer_Indice_Al_Final(tensor_estres)[1][1]);
		serie_espacial[25] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[26] = Util.Vector_Componentes_Iguales(N,0.0);
		serie_espacial[27] = Util.Vector_Componentes_Iguales(N,0.0);
		serie_espacial[28] = Util.Vector_Componentes_Iguales(N,0.0);	  
		serie_espacial[29] = Util.Vector_Componentes_Iguales(N,0.0);
		serie_espacial[30] = Util.Vector_Componentes_Iguales(N,0.0);
		serie_espacial[31] = Util.Vector_Componentes_Iguales(N,0.0);	  
		serie_espacial[32] = Util.Vector_Componentes_Iguales(N,0.0);	  
		serie_espacial[33] = Util.Copiar_Vector(Util.Transpuesta(estres_prin)[0]);		  
		serie_espacial[34] = Util.Copiar_Vector(Util.Transpuesta(estres_prin)[1]);
		serie_espacial[35] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[36] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[37] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[38] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[39] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[40] = Util.Vector_Componentes_Iguales(N,0.0);	
		serie_espacial[42] = Util.Copiar_Vector(h);
		
		  
		return serie_espacial;
	}
	

	
		public String[] Nombre_Serie_Espacial()
	{
		String[] nombres = new String[numero_variables_guardar];
		
		nombres[0]  = "posicion.en.x";			
		nombres[1]  = "posicion.en.y";
		nombres[2]  = "elongacion.en.x";
		nombres[3]  = "elongacion.en.y";                       
		nombres[4]  = "posicion.particulas.con.respecto.centro.de.masa.x";         
		nombres[5]  = "posicion.particulas.con.respecto.centro.de.masa.y";   	
		nombres[6]  = "velocidad.en.x";			
		nombres[7]  = "velocidad.en.y";       
		nombres[8]  = "velocidad.particulas.con.respecto.centro.de.masa.x";       
		nombres[9]  = "velocidad.particulas.con.respecto.centro.de.masa.y";         
		nombres[10] = "aceleracion.en.x";            
		nombres[11] = "aceleracion.en.y"; 
		nombres[12] = "densidad";		
		nombres[13] = "tasa.densidad";			
		nombres[14] = "energia_elastica .helmontz";
		nombres[15] = "tasa.energia_elastica .helmontz";
		nombres[16] = "energia_elastica .interna";
		nombres[17] = "tasa.energia_elastica .interna";
		nombres[18] = "energia_elastica .cinetica";
		nombres[19] = "deformacion.x.x";            
		nombres[20] = "deformacion.x.y";            
		nombres[21] = "deformacion.y.y";          
		nombres[22] = "estres.x.x";
		nombres[23] = "estres.x.y";            
		nombres[24] = "estres.y.y"; 
		nombres[25] = "presion";		
		nombres[26] = "deviatorico.x.x";			
		nombres[27] = "deviatorico.x.y";
		nombres[28] = "deviatorico.y.y";
		nombres[29] = "tasa.deviatorico.x.x";	            
		nombres[30] = "tasa.deviatorico.x.y";            
		nombres[31] = "tasa.deviatorico.y.y";          
		nombres[32] = "angulo.estres.principal";
		nombres[33] = "estres.principal.primero";            
		nombres[34] = "estres.principal.segundo";     
		nombres[35] = "suma.kernel";		
		nombres[36] = "suma.kernel.con.consistencia";			
		nombres[37] = "suma.gradiente.kernel.x";
		nombres[38] = "suma.gradiente.kernel.y";
		nombres[39] = "suma.gradiente.kernel.x.con.consistencia";           
		nombres[40] = "suma.gradiente.kernel.y.con.consistencia"; 
		nombres[41] = "masa.particula";           
		nombres[42] = "distancia.suavizado";   

		return nombres;
	}
}
