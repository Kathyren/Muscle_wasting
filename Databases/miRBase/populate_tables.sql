use mirbase;
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/confidence_score.txt' INTO TABLE confidence_score ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/dead_mirna.txt' INTO TABLE dead_mirna ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/literature_references.txt' INTO TABLE literature_references ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mature.txt' INTO TABLE mature ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mature_database_links.txt' INTO TABLE mature_database_links ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mature_database_url.txt' INTO TABLE mature_database_url ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna.txt' INTO TABLE mirna ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_2_prefam.txt' INTO TABLE mirna_2_prefam ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_chromosome_build.txt' INTO TABLE mirna_chromosome_build ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_context.txt' INTO TABLE mirna_context ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_database_links.txt' INTO TABLE mirna_database_links ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_database_url.txt' INTO TABLE mirna_database_url ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_literature_references.txt' INTO TABLE mirna_literature_references ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_mature.txt' INTO TABLE mirna_mature ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_pre_mature.txt' INTO TABLE mirna_pre_mature ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_prefam.txt' INTO TABLE mirna_prefam ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/mirna_species.txt' INTO TABLE mirna_species ;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRBase/organisms.txt' INTO TABLE organisms ;

-- Parsed
use mirnadbs;
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/parsed/mirTarBase_miRTarBase_Target.txt' INTO TABLE binding;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/parsed/mirWalk_dre_miRWalk_3UTR.txt' INTO TABLE binding;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/parsed/mirWalk_hsa_miRWalk_5UTR.txt' INTO TABLE binding;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/parsed/mirWalk_mmu_miRWalk_3UTR.txt' INTO TABLE binding;
use mirnadbs;
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/parsed/miRDB_miRDB_v6.0_prediction_result.txt' INTO TABLE binding;
use mirnadbs;
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/parsed/genes.txt' INTO TABLE gene_bank;



use mirnadbs;
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/targetscan/seeds.txt' INTO TABLE mirna_seeds ;

use mirnadbs;
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/miRNATissueAtlas2/miRNATissueAtlas2.txt' INTO TABLE mirna_tissues ;

-- Adding experiments for Home sample (dunhill)
use mirnadbs;
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/dunhill/experiments.txt' INTO TABLE experiments (experiment_id,source,description,author,link ) ;
-- Adding the unique organs
SET GLOBAL local_infile=1;
LOAD DATA LOCAL INFILE 'C:/Users/crtuser/Documents/PhD/Project/Databases/organs.txt' INTO TABLE organs (organ) ;


-- Adding

