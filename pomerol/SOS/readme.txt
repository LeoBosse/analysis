Version 5.0 : 01/06/2010
-----------

- Complete reprise du parametrage : utilisation de "-keyword value".

- Introduction d'une fonction SOS_ANGLES.F qui gere la definition des angles de Gauss
  et l'utilisation d'angles utilisateurs. Elle gere également les ordres limites des developpements.
	
- Definition de tableaux (luminances, fonctions de phase) sur des domaines plus importants 
  que le domaine utilisé : permet un usage du code avec une possibilité de modification du 
  nombre d'angles utilisés sans recompilation (et l'ajout d'angles utilisateurs variables).
  
- Adaptation de la routine SOS_AEROSOLS.F pour permettre l'utilisation d'un fichier 
  de fonctions de phase externes (utile pour des particules non-spheriques).
  
- Adaptation de la routine SOS.F pour permettre le calcul des transmissions diffuses.
	    
- Tracabilité des constantes par copie du fichier SOS.h sur l'espace de compilation 
  et sur l'espace des resultats (modification du main_SOS.ksh).

- Mise a jour du Manuel Utilisateur + version anglaise.

       
######################################################       
       
Version 4.1 : 08/09/2008
-----------

- Intègre une prise en compte des modeles bimodaux (Dubovik) avec correctifs / version 4.0.

- Evolutions depuis la version 4.0 :
     - Adaptation de la routine SOS_AEROSOLS.F (VERSION:2.1):
        *  Correction de depassements de zone d'indentation pour compilation f77.
	
        * Suppression de la sortie des calculs pour fin de fichier de MIE atteinte.
          Car si c'est le cas, il s'agit d'une erreur : le fichier de MIE est 
          incomplet (probable manque d'espace disque à sa génération).
          --> Gestion du cas d'erreur.
	  
        * Correction pour le calcul des proportions de composants des modèles bimodaux
          Dubovik :
	  Cas du calcul à partir des rapports d'épaisseur optique des deux modes
          (fin et grossier) et de l'épaisseur optique totale pour une longueur
           d'onde de référence 
	  --> utilisation de l'indice de réfraction des particules pour la longueur 
	  d'onde de référence (uniquement).
	    
     - Evolution du script main_SOS.ksh pour l'appel de la fonction SOS_AEROSOLS.exe
       (gestion des arguments)


