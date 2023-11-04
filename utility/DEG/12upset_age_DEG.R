library(UpSetR)
input=c('BC&Cone&Rod'=26,'AC'=583,'BC&HC&RGC'=9,'AC&BC&MG'=17,'AC&BC&Cone&HC&MG'=3,'BC&HC'=200,'AC&BC&HC&MG&RGC'=8,'Cone'=274,'BC'=704,'AC&HC&MG'=10,'AC&BC&HC'=181,'AC&BC&Cone&MG&Rod'=9,'AC&Rod'=129,'AC&BC&RGC'=19,'AC&Cone&Rod'=22,'MG&RGC'=4,'AC&BC&Cone&HC&Rod'=32,'BC&Cone'=29,'Cone&Rod'=145,'BC&HC&Rod'=63,'HC&MG&RGC'=1,'BC&Cone&MG&RGC'=1,'Cone&HC&RGC'=3,'AC&BC&HC&Rod'=89,'AC&BC&Rod'=113,'AC&RGC'=31,'HC&Rod'=99,'AC&Cone&RGC&Rod'=6,'HC&RGC'=27,'BC&Cone&MG&RGC&Rod'=1,'HC&MG'=23,'AC&BC&HC&MG'=31,'AC&BC&Cone&HC&MG&RGC&Rod'=7,'Cone&HC'=36,'BC&Cone&RGC'=1,'BC&RGC&Rod'=8,'MG&Rod'=28,'HC&MG&Rod'=2,'AC&BC&Cone&Rod'=40,'Cone&MG&RGC'=1,'BC&Cone&HC&MG&Rod'=1,'BC&HC&RGC&Rod'=3,'AC&Cone&MG'=1,'AC&BC&MG&RGC'=3,'BC&HC&MG'=13,'AC&BC&MG&Rod'=18,'Cone&HC&Rod'=12,'MG'=150,'Cone&MG&Rod'=3,'AC&BC&Cone&MG'=1,'AC&Cone&HC&MG&RGC'=1,'AC&BC&Cone&MG&RGC&Rod'=3,'BC&MG&RGC&Rod'=1,'AC&Cone&HC&MG&Rod'=3,'BC&Rod'=207,'AC&BC'=232,'AC&Cone&MG&Rod'=1,'AC&Cone&HC'=8,'AC&BC&Cone&HC&RGC'=1,'AC&Cone&RGC'=3,'AC&BC&Cone&HC&RGC&Rod'=15,'AC&Cone&HC&Rod'=10,'AC&HC&MG&Rod'=4,'Cone&MG&RGC&Rod'=3,'AC&BC&Cone&HC'=17,'BC&MG'=35,'AC&BC&Cone&HC&MG&RGC'=2,'Rod'=1491,'Cone&HC&RGC&Rod'=3,'BC&MG&RGC'=1,'AC&HC&RGC'=10,'BC&Cone&HC&Rod'=9,'BC&RGC'=27,'AC&MG&Rod'=7,'BC&MG&Rod'=15,'BC&Cone&HC'=11,'AC&BC&HC&RGC&Rod'=16,'Cone&MG'=3,'Cone&RGC'=13,'AC&BC&HC&MG&RGC&Rod'=13,'AC&Cone'=20,'AC&MG'=11,'BC&Cone&MG&Rod'=2,'AC&BC&MG&RGC&Rod'=6,'AC&BC&HC&RGC'=23,'AC&BC&HC&MG&Rod'=32,'RGC'=344,'HC&MG&RGC&Rod'=1,'HC&RGC&Rod'=5,'AC&BC&Cone&HC&MG&Rod'=7,'AC&BC&Cone'=21,'AC&BC&Cone&RGC&Rod'=6,'RGC&Rod'=70,'AC&HC'=136,'BC&Cone&RGC&Rod'=4,'Cone&RGC&Rod'=10,'HC'=960,'AC&HC&Rod'=36,'BC&Cone&MG'=2,'BC&HC&MG&Rod'=8,'AC&BC&RGC&Rod'=18,'AC&RGC&Rod'=11,'AC&HC&RGC&Rod'=2,'AC&Cone&HC&RGC&Rod'=1)
pdf("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/upset_list_all_age_all.pdf",width=10)
upset(fromExpression(input), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = "DEG per cell type",
      sets.bar.color="blue"
      )
dev.off()

input=c('BC&Cone&Rod'=13,'AC'=445,'BC&HC&RGC'=2,'AC&BC&MG'=8,'BC&HC'=78,'AC&BC&HC&MG&RGC'=3,'Cone'=163,'BC'=343,'AC&HC&MG'=6,'AC&BC&HC'=56,'AC&BC&Cone&MG&Rod'=5,'AC&Rod'=82,'AC&BC&RGC'=8,'AC&Cone&Rod'=16,'MG&RGC'=1,'AC&BC&Cone&HC&Rod'=18,'BC&Cone'=16,'Cone&Rod'=75,'BC&HC&Rod'=13,'HC&MG&RGC'=1,'BC&Cone&MG&RGC'=1,'AC&BC&HC&Rod'=28,'AC&BC&Rod'=57,'AC&RGC'=18,'HC&Rod'=28,'HC&RGC'=1,'AC&Cone&RGC&Rod'=5,'HC&MG'=14,'AC&BC&HC&MG'=7,'AC&BC&Cone&HC&MG&RGC&Rod'=6,'Cone&HC'=8,'BC&Cone&RGC'=1,'BC&RGC&Rod'=1,'MG&Rod'=17,'HC&MG&Rod'=2,'AC&BC&Cone&Rod'=35,'Cone&MG&RGC'=1,'BC&Cone&HC&MG&Rod'=1,'AC&Cone&MG'=1,'AC&BC&MG&RGC'=1,'BC&HC&MG'=4,'AC&BC&MG&Rod'=8,'Cone&HC&Rod'=3,'MG'=96,'Cone&MG&Rod'=3,'AC&Cone&HC&MG&RGC'=1,'AC&BC&Cone&MG&RGC&Rod'=3,'BC&MG&RGC&Rod'=1,'AC&Cone&HC&MG&Rod'=2,'BC&Rod'=108,'AC&BC'=135,'AC&Cone&MG&Rod'=1,'AC&Cone&HC'=5,'AC&Cone&RGC'=2,'AC&BC&Cone&HC&RGC&Rod'=5,'AC&Cone&HC&Rod'=7,'AC&HC&MG&Rod'=1,'Cone&MG&RGC&Rod'=2,'AC&BC&Cone&HC'=7,'BC&MG'=14,'AC&BC&Cone&HC&MG&RGC'=1,'Rod'=756,'Cone&HC&RGC&Rod'=2,'BC&MG&RGC'=1,'AC&HC&RGC'=5,'BC&RGC'=12,'BC&Cone&HC&Rod'=4,'AC&MG&Rod'=5,'AC&BC&HC&RGC&Rod'=1,'BC&Cone&HC'=1,'BC&MG&Rod'=9,'Cone&MG'=2,'Cone&RGC'=3,'AC&BC&HC&MG&RGC&Rod'=3,'AC&Cone'=9,'AC&MG'=5,'AC&BC&MG&RGC&Rod'=1,'AC&BC&HC&RGC'=5,'BC&Cone&MG&Rod'=2,'AC&BC&HC&MG&Rod'=5,'RGC'=112,'HC&RGC&Rod'=1,'AC&BC&Cone&HC&MG&Rod'=5,'AC&BC&Cone'=12,'AC&BC&Cone&RGC&Rod'=3,'RGC&Rod'=10,'AC&HC'=75,'Cone&RGC&Rod'=4,'HC'=469,'AC&HC&Rod'=14,'BC&Cone&MG'=1,'AC&RGC&Rod'=8,'AC&BC&RGC&Rod'=9,'AC&HC&RGC&Rod'=2)

#pdf("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch/upset_list_up_rd_new.pdf",width=15)
pdf("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/upset_list_all_age_up.pdf",width=10)
upset(fromExpression(input),
      nintersects = 5000,
      nsets = 20,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000,
      text.scale = c(2, 1.5, 2, 1.15, 2, 0),
      point.size = 3.5,
      line.size = 1
      mainbar.y.label = "DEG intersections",
      sets.x.label = "DEG per cell type",
      sets.bar.color="blue"

      )
dev.off()

input=c('BC&Cone&Rod'=13,'AC'=157,'BC&HC&RGC'=7,'AC&BC&MG'=9,'AC&BC&Cone&HC&MG'=3,'BC&HC'=125,'AC&BC&HC&MG&RGC'=5,'BC'=371,'Cone'=133,'AC&HC&MG'=4,'AC&BC&HC'=125,'AC&BC&Cone&MG&Rod'=4,'AC&Rod'=42,'AC&BC&RGC'=11,'AC&Cone&Rod'=6,'MG&RGC'=2,'BC&Cone'=13,'AC&BC&Cone&HC&Rod'=14,'Cone&Rod'=68,'BC&HC&Rod'=47,'Cone&HC&RGC'=2,'AC&BC&HC&Rod'=61,'AC&BC&Rod'=56,'AC&RGC'=12,'HC&Rod'=40,'AC&Cone&RGC&Rod'=1,'HC&RGC'=19,'BC&Cone&MG&RGC&Rod'=1,'HC&MG'=10,'AC&BC&HC&MG'=24,'AC&BC&Cone&HC&MG&RGC&Rod'=1,'Cone&HC'=15,'BC&RGC&Rod'=4,'MG&Rod'=9,'AC&BC&Cone&Rod'=5,'BC&HC&RGC&Rod'=2,'AC&BC&MG&RGC'=2,'BC&HC&MG'=9,'AC&BC&MG&Rod'=10,'Cone&HC&Rod'=8,'MG'=60,'AC&BC&Cone&MG'=1,'AC&Cone&HC&MG&Rod'=1,'BC&Rod'=98,'AC&BC'=96,'AC&Cone&HC'=3,'AC&BC&Cone&HC&RGC'=1,'AC&Cone&RGC'=1,'AC&BC&Cone&HC&RGC&Rod'=10,'AC&Cone&HC&Rod'=3,'AC&HC&MG&Rod'=3,'Cone&MG&RGC&Rod'=1,'AC&BC&Cone&HC'=10,'BC&MG'=20,'AC&BC&Cone&HC&MG&RGC'=1,'Rod'=800,'Cone&HC&RGC&Rod'=1,'AC&HC&RGC'=5,'BC&Cone&HC&Rod'=5,'BC&RGC'=12,'AC&MG&Rod'=2,'BC&MG&Rod'=6,'AC&BC&HC&RGC&Rod'=15,'BC&Cone&HC'=10,'Cone&MG'=1,'Cone&RGC'=9,'AC&BC&HC&MG&RGC&Rod'=10,'AC&MG'=4,'AC&Cone'=7,'AC&BC&MG&RGC&Rod'=5,'AC&BC&HC&RGC'=18,'AC&BC&HC&MG&Rod'=27,'RGC'=262,'HC&RGC&Rod'=1,'AC&BC&Cone&HC&MG&Rod'=2,'AC&BC&Cone'=9,'AC&BC&Cone&RGC&Rod'=3,'RGC&Rod'=52,'AC&HC'=59,'BC&Cone&RGC&Rod'=4,'Cone&RGC&Rod'=6,'HC'=554,'AC&HC&Rod'=18,'BC&Cone&MG'=1,'BC&HC&MG&Rod'=8,'AC&BC&RGC&Rod'=9,'AC&RGC&Rod'=3,'AC&Cone&HC&RGC&Rod'=1)

#pdf("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch/upset_list_down_rd_new.pdf",width=15)
pdf("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/upset_list_all_age_down.pdf",width=10)

upset(fromExpression(input),
      nintersects = 5000,
      nsets = 20,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000,
      text.scale = c(2, 1.5, 2, 1.15, 2, 0),
      point.size = 3.5,
      line.size = 1
      mainbar.y.label = "DEG intersections",
     sets.x.label = "DEG per cell type",
	sets.bar.color="blue" 
      )
dev.off()
