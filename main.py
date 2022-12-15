# -*- coding: utf-8 -*-
"""
Projeto nº2 - Fotogrametria Digital - Mestrado Eng. Geoespacial 2022-23
Fotogrametria analítica - cálculo de orientação externa e coordenadas objeto de
pontos novos num modelo esteroescópico.

@author: Mário Amaral
"""
import numpy as np

def main():
    """O programa escreve os valores da orientação externa da foto 0006 e das
        coordenadas dos novos pontos para um ficheiro de resultados em txt na
        mesma diretoria em que o programa é invocado.
    """
    
    output_filename = 'resultados.txt' # nome do ficheiro para escrita dos resultados finais (na mesma diretoria do programa)
    
    # Dados das coordenadas foto, coordenadas objeto e orientação das fotos
    # são guardados em dicionários python para fácil referência.
    
    coord_photo = {'0006' : {'1000' : np.array([[-1.026], 
                                                [68.670]]),
                             '2000' : np.array([[44.310], 
                                                [65.634]]),
                             '3000' : np.array([[2.526], 
                                                [-70.734]]),
                             '4000' : np.array([[41.166], 
                                                [-61.722]]),
                             '5000' : np.array([[-1.962], 
                                                [7.386]]),	
                             '6000' : np.array([[36.534], 
                                                [0.642]]),
                             'A' : np.array([[38.898], 
                                             [-34.266]]),
                             'B' : np.array([[34.268], 
                                             [-47.658]]),
                             'C' : np.array([[41.986], 
                                             [-50.164]])
                              },
                    
                    '0007' : {'A' : np.array([[4.866], 
                                              [-38.166]]),
                              'B' : np.array([[0.606], 
                                              [-51.966]]),
                              'C' : np.array([[8.490], 
                                              [-54.402]])
                              }
                    }
    
    coord_obj = {'1000' : np.array([[-90158.653], 
                                    [-100565.532], 
                                    [93.799]]),
                 '2000' : np.array([[-89794.832], 
                                    [-100478.941], 
                                    [89.069]]),
                 '3000' : np.array([[-89788.239], 
                                    [-101679.518], 
                                    [81.974]]),
                 '4000' : np.array([[-89514.366], 
                                    [-101489.343], 
                                    [123.038]]),
                 '5000' : np.array([[-90019.770], 
                                    [-101048.356], 
                                    [93.223]]),
                 '6000' : np.array([[-89696.832], 
                                    [-101008.761], 
                                    [95.277]]),
                 'A' : np.zeros((3,1)), #desconhecidos à partida
                 'B' : np.zeros((3,1)), #desconhecidos à partida
                 'C' : np.zeros((3,1)), #desconhecidos à partida
                 }
    
    orient_photo = {'0006' : np.zeros((6,1)), #desconhecida à partida
                    '0007' : np.array([[-89711.171], 
                                       [-100989.573], 
                                       [1094.305], 
                                       [0.1208], 
                                       [0.2895], 
                                       [15.7612]])    
                    }
    c = 120 # distância focal (em mm)
    
    ###
    ### 0. Funções auxiliares para cálculo de coordenadas com base nas ECOL
    ###
    
    def get_R(OE):
        """ Função retorna a matriz de rotação com base na orientação externa (OE) 
            de uma foto.  
            R[0,0] = r11, ..., R[3,3] = r33
            Ângulos convertidos para radianos
        """
        
        #ângulos convertidos para radianos
        omega = np.radians(OE[3,0]) 
        phi = np.radians(OE[4,0]) 
        kappa = np.radians(OE[5,0]) 
        
        return np.array([[np.cos(phi) * np.cos(kappa), - np.cos(phi) * np.sin(kappa), np.sin(phi)],
                         [np.cos(omega) * np.sin(kappa) + np.sin(omega) * np.sin(phi) * np.cos(kappa), np.cos(omega) * np.cos(kappa) - np.sin(omega) * np.sin(phi) * np.sin(kappa), - np.sin(omega) * np.cos(phi)],
                         [np.sin(omega) * np.sin(kappa) - np.cos(omega) * np.sin(phi) * np.cos(kappa), np.sin(omega) * np.cos(kappa) + np.cos(omega) * np.sin(phi) * np.sin(kappa), np.cos(omega) * np.cos(phi)]])
    
    
    def get_N(P, OE, R):
        """ Função retorna vector do numerador das ECOL calculado para um determinado 
            ponto com coordenadas objeto P, orientação externa OE e matriz de 
            rotação R.
            Nx = N[0], Ny = N[1]
        """
        X,Y,Z = P[:,0] 
        XO,YO,ZO = OE[:3,0]
        return np.array([R[0,0]*(X-XO) + R[1,0]*(Y-YO) + R[2,0]*(Z-ZO),
                         R[0,1]*(X-XO) + R[1,1]*(Y-YO) + R[2,1]*(Z-ZO)])
    
    
    def get_D(P, OE, R):
        """ Função retorna denominadar das ECOL calculado para um determinado 
            ponto com coordenadas objeto P, orientação externa OE e matriz de 
            rotação R.
        """
        X,Y,Z = P[:,0]
        XO,YO,ZO = OE[:3,0]
        return R[0,2]*(X-XO) + R[1,2]*(Y-YO) + R[2,2]*(Z-ZO)
                   
    
    def get_coord_photo(OE, P):
        """ Função calcula as coordenadas foto para um determinado ponto com coord
            objeto P, com base na orientação externa OE.
        """
        R = get_R(OE)
        Nx, Ny = get_N(P, OE, R)
        D = get_D(P, OE, R)
        
        x_o = -c * (Nx/D) 
        y_o = -c * (Ny/D) 
        
        return np.array([[x_o],
                         [y_o]])
    
    
    def print_coord(output_filename):
        """ Escreve para ficheiro as coordenadas da orientação externa e coord objeto dos pontos novos com a designação
            dada por label
        """
        photo_id = '0006'
        list_point_ids = ['A','B','C']
            
        with open(output_filename, 'w') as file:
            
            XO, YO, ZO, omega, phi, kappa = orient_photo[photo_id][:,0]
            output = 'Orientação externa da foto ' + photo_id + \
                  ': X= {0:1.2f}, Y= {1:1.2f}, Z= {2:1.2f}, omega = {3:1.2f}, phi = {4:1.2f}, kappa = {5:1.2f}'.format(XO,YO,ZO,omega,phi,kappa)
            
            file.write(output + "\n")
            
            for point_id in list_point_ids:
                X, Y, Z = coord_obj[point_id][:,0]
                output = 'Coordenadas obj do ponto ' + point_id + ': X= {0:1.2f}, Y= {1:1.2f}, Z= {2:1.2f} '.format(X,Y,Z)
                file.write(output + "\n")
    
    ###
    ### 1. Determinação da orientação externa da foto 0006
    ###
    
    def get_delta_l_orient_ext(OE, photo_id, point_id):
        """ Função calcula o vector delta_l das observações para um determinado
            valor de orientação externa OE, com vista à determinação da OE pelo
            método dos mínimos quadrados.
            p_o: coordenadas foto x_o,y_o dadas pelas ECOL para o ponto identificado 
                 pelo point_id, com base na OE fornecida
            p_barr: coordenadas x_barr, y_barrfoto medidas e fornecidas como dados 
                 para o projeto
            photo_id: identificação da foto a que dizem respeito as coordenadas
        """        
        p_o = get_coord_photo(OE, coord_obj[point_id])
        p_barr = coord_photo[photo_id][point_id]                
            
        return p_o - p_barr
        
    
    def get_A_orient_ext(OE, point_id):
        """ Função calcula o Jacobiano A para determinação da orientação externa da
            foto pelo método dos mínimos quadrados. Os valores são determinados para
            uma dada orientação externa OE e para um dado ponto de coordenadas obj
            conhecidas (identificado através do point_id).
            Retorna um vector de 2 linhas para x, y e de 6 colunas para as 6 incógnitas
            XO, YO, ZO, omega, phi, kappa. Um par de linhas por ponto.
        """       
        P = coord_obj[point_id]
        R = get_R(OE)
        Nx,Ny = get_N(P,OE,R)
        D = get_D(P, OE, R)
        
        X,Y,Z = P[:,0]    
        XO,YO,ZO = OE[:3,0]
        
        kappa = np.radians(OE[5,0])
    
        dx_dXO = -(c/D**2) * (R[0,2]*Nx - R[0,0]*D)
        dy_dXO = -(c/D**2) * (R[0,2]*Ny - R[0,1]*D)
        
        dx_dYO = -(c/D**2) * (R[1,2]*Nx - R[1,0]*D)
        dy_dYO = -(c/D**2) * (R[1,2]*Ny - R[1,1]*D)
        
        dx_dZO = -(c/D**2) * (R[2,2]*Nx - R[2,0]*D)        
        dy_dZO = -(c/D**2) * (R[2,2]*Ny - R[2,1]*D)
        
        dx_domega = -(c/D) * (((Y-YO)*R[2,2] - (Z-ZO)*R[1,2])*(Nx/D) - (Y-YO)*R[2,0] + (Z-ZO)*R[1,0])
        dy_domega = -(c/D) * (((Y-YO)*R[2,2] - (Z-ZO)*R[1,2])*(Ny/D) - (Y-YO)*R[2,1] + (Z-ZO)*R[1,1])
        
        dx_dphi = (c/D) * ((Nx*np.cos(kappa) - Ny*np.sin(kappa))*(Nx/D) + D*np.cos(kappa))
        dy_dphi = (c/D) * ((Nx*np.cos(kappa) - Ny*np.sin(kappa))*(Ny/D) - D*np.sin(kappa))
        
        dx_dkappa = -(c/D) * Ny
        dy_dkappa =  (c/D) * Nx
        
        return np.array([[dx_dXO, dx_dYO, dx_dZO, dx_domega, dx_dphi, dx_dkappa],
                          [dy_dXO, dy_dYO, dy_dZO, dy_domega, dy_dphi, dy_dkappa]])
    
    
    ## Valores aproximados iniciais para a orientação externa da foto 0006 são calculados
    ## a partir das coordenadas dos PFs medidos no terro. A coordenada ZO assume-se com
    ## valor inicial igual à da foto 0007. Omega, phi assumem-se nulos e o ângulo kappa
    ## assume-se igual ao ângulo kappa da fotografia 0007
    
    
    X_avg = np.average([coord_obj['1000'][0,0],
                       coord_obj['2000'][0,0],
                       coord_obj['3000'][0,0],
                       coord_obj['4000'][0,0],
                       coord_obj['5000'][0,0],
                       coord_obj['6000'][0,0]])
    
    Y_avg = np.average([coord_obj['1000'][1,0],
                       coord_obj['2000'][1,0],
                       coord_obj['3000'][1,0],
                       coord_obj['4000'][1,0],
                       coord_obj['5000'][1,0],
                       coord_obj['6000'][1,0]])
    
    OE = np.array([[X_avg],
                   [Y_avg], 
                   [1094.], 
                   [0.], 
                   [0.], 
                   [15.]])
    
    photo_id = '0006'
    point1 = '1000'
    point2 = '4000'
    point3 = '6000'
    
    delta_X_track = np.zeros((6,1)) #array para registar iterações de delta_X
    OE_track = np.zeros((6,1)) #array para registar iterações de OE aproximado
    
    for i in range(100): # definimos arbitrariamente um limite de 100 iterações, sabendo que há convergência antes disso
    
        A = np.vstack((get_A_orient_ext(OE, point1),
                       get_A_orient_ext(OE, point2),
                       get_A_orient_ext(OE, point3)))
        
        delta_l = np.vstack((get_delta_l_orient_ext(OE, photo_id, point1),
                             get_delta_l_orient_ext(OE, photo_id, point2), 
                             get_delta_l_orient_ext(OE, photo_id, point3)))
        
        N = np.transpose(A) @ A
        h = np.transpose(A) @ delta_l
        
        delta_X = np.linalg.inv(N) @ h 
        
        if abs(delta_X[0,0]) < 0.1: break # o for loop é interrompido para um |delta_X| inferior a 10cm
        
        OE = OE - delta_X
          
        delta_X_track = np.hstack((delta_X_track, delta_X))
        OE_track = np.hstack((OE_track, OE))
        
    orient_photo['0006'] = OE ## variável guardada para output no fim do programa
    
    ## Os valores de delta_X não convergem para zero e infelizmente são diferentes
    ## conforme diferentes combinações de pontos.
    ## Não consigo entender onde poderá estar o problema, depois de múltiplos testes
    ## e revisões do código. 
    ## Como não consegui obter convergência nos valores de OE para a foto 0006, para a
    ## etapa seguinte, optei por utilizar os valores de orientação externa da foto
    ## 0006 obtidos no exercício com o Photomod.
    
        
    ###
    ### 2. Determinação das coordenadas objeto dos novos pontos A, B e C por interseção espacial
    ### direta
    ###
    
    ###
    ### 2.1 Funções auxiliares para cáclulo das coordenadas com base nas ECOL e
    ###     determinação das coordenadas aproximadas iniciais
    ###
    
    def get_K(OE):
        """ Função auxiliar para cálculo dos coeficientes Kx e Ky necessários para 
            determinar as coordenadas objeto iniciais para um novo ponto. Estes 
            coeficientes são calculados com base na orientação externa de uma das fotos (OE)
        """
        R = get_R(OE)
        
        x = OE[0,0]
        y = OE[1,0]
        
        Kx = (R[0,0]*x + R[0,1]*y - R[0,2]*c)/(R[2,0]*x + R[2,1]*y - R[2,2]*c)
        Ky = (R[1,0]*x + R[1,1]*y - R[1,2]*c)/(R[2,0]*x + R[2,1]*y - R[2,2]*c)
    
        return np.array([Kx, Ky])
    
    
    def get_P_init(photo_key1, photo_key2):
        """ Função calcula as coordenadas objeto aproximadas iniciais tendo por base
            as orientações externas das duas fotos que constituem o modelo onde o
            ponto novo é localizado.
        """ 
        
        OE1 = orient_photo[photo_key1]
        OE2 = orient_photo[photo_key2]
        
        XO1, YO1, ZO1 = OE1[:3,0]
        XO2, YO2, ZO2 = OE2[:3,0]
        
        K1x, K1y = get_K(OE1)
        K2x, K2y = get_K(OE2)
        
        Z = (XO2 - ZO2*K2x + ZO1*K1x - XO1) / (K1x - K2x)
        X = ((XO1 + (Z - ZO1)*K1x) + ((XO2 + (Z - ZO2)*K2x)))/2.
        Y = ((YO1 + (Z - ZO1)*K1y) + ((YO2 + (Z - ZO2)*K2y)))/2.
        
        return np.array([[X],
                         [Y],
                         [Z]])        
    
    
    def get_delta_l_new_P(photo_id, new_point_id, P):
        """ Função calcula o vector delta_l das observações para um determinado
            valor de coordenadas objeto aproximadas do novo ponto (P), com vista à
            sua determinação pelo método dos mínimos quadrados.
            p_o: coordenadas foto x_o,y_o dadas pelas ECOL para o valor aproximado
                 do ponto (P)
            p_barr: coordenadas x_barr, y_barrfoto medidas na foto, fornecidas 
                    como dados para o projeto 
            photo_id: identificação da foto a que dizem respeito as coordenadas
        """        
        p_o = get_coord_photo(orient_photo[photo_id], P)
        p_barr = coord_photo[photo_id][new_point_id]                
            
        return p_o - p_barr
    
    
    def get_A_new_P(photo_id, P):
        """ Função calcula o Jacobiano A para determinação das coordenadas objeto do
            novo ponto pelo método dos mínimos quadrados. Os valores são determinados 
            para o ponto aproximado P, tendo por base as orientações externas das
            fotos no modelo esteroscópico.
            Retorna um vector de 2 linhas para x, y e de 3 colunas para as 3 incógnitas
            X, Y, Z. Um par de linhas por foto.
        """  
        OE = orient_photo[photo_id]
        
        R = get_R(OE)
        Nx, Ny = get_N(P, OE, R)
        D = get_D(P, OE, R)
    
        dx_dX = -(c/D**2) * (R[0,0]*D - R[0,2]*Nx)
        dx_dY = -(c/D**2) * (R[1,0]*D - R[1,2]*Nx)
        dx_dZ = -(c/D**2) * (R[2,0]*D - R[2,2]*Nx)        
                
        dy_dX = -(c/D**2) * (R[0,1]*D - R[0,2]*Ny)
        dy_dY = -(c/D**2) * (R[1,1]*D - R[1,2]*Ny)
        dy_dZ = -(c/D**2) * (R[2,1]*D - R[2,2]*Ny)        
        
        return np.array([[dx_dX, dx_dY, dx_dZ],
                         [dy_dX, dy_dY, dy_dZ]])
    
    ## 2.2 Determinação das coordenadas objeto do novo ponto A
    
    ## Atendendo a que não se consegue convergência na determinação da orientação externa da 
    ## foto 0006 no passo anterior, utiliza-se aqui o valor de OE obtido através do Photomod.
    
    ## Em cada caso a convergência obteve-se rapidamente, dentro de 20 iterações
    
    P_aprox = get_P_init('0006', '0007')
    
    delta_X_track = np.zeros(P_aprox.shape) #array para registar iterações de delta_X
    P_aprox_track = np.zeros(P_aprox.shape) #array para registar iterações das coordenadas do ponto aproximado
    
    # list_new_points = [coord_obj['A'], coord_obj['B'], coord_obj['C']]
    
    list_new_point_ids = ['A', 'B', 'C']
    
    for point_id in list_new_point_ids:
        
        for i in range(100): # definimos arbitrariamente um limite de 100 iterações, sabendo que há convergência antes disso
        
            A = np.vstack((get_A_new_P('0006', P_aprox),
                           get_A_new_P('0007', P_aprox)))
            
            delta_l = np.vstack((get_delta_l_new_P('0006', point_id, P_aprox),
                                 get_delta_l_new_P('0007', point_id, P_aprox)))
            
            N = np.transpose(A) @ A
            h = np.transpose(A) @ delta_l
            
            delta_X = np.linalg.inv(N) @ h 
            
            if abs(delta_X[0,0]) < 0.1: break # o for loop é interrompido para um |delta_X| inferior a 10cm
            
            P_aprox = P_aprox - delta_X
            
            delta_X_track = np.hstack((delta_X_track, delta_X))
            P_aprox_track = np.hstack((P_aprox_track, P_aprox))
        
        coord_obj[point_id]= P_aprox
    
    print_coord(output_filename)


if __name__ == "__main__":
    main()