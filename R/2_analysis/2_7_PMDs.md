---
title: "2_7_PMDs"
author: "Victor Yuan"
date: "23/09/2019"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: false
    theme: spacelab
editor_options: 
  chunk_output_type: console
---


# Setup

## Libraries



```r
# libraries and data
library(tidyverse); theme_set(theme_bw())
library(DT)
```

## Data


```r
pDat <- readRDS('../../data/main/interim/2_3_pDat_contam.rds')
pDat <- pDat %>%
  mutate(Tissue = case_when(
    !(Tissue %in% c('Villi', 'Villi maternal', 'Syncytiotrophoblast')) ~ paste(Tissue, 'cs'),
    Tissue == 'Syncytiotrophoblast' ~ 'Trophoblasts enz',
    TRUE ~ Tissue
  )) 

# raw methylation data
betas <- readRDS('../../data/main/interim/1_4_betas_noob_filt.rds')

mset_noob <- readRDS('../../data/main/interim/1_4_mset_noob.rds') # for mvals
colnames(mset_noob) <- colnames(betas) <- pDat$Sample_Name
```

```
## Loading required package: minfi
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unsplit, which,
##     which.max, which.min
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## The following object is masked from 'package:purrr':
## 
##     reduce
```

```
## The following object is masked from 'package:grDevices':
## 
##     windows
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## The following object is masked from 'package:dplyr':
## 
##     count
```

```
## Loading required package: BiocParallel
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following object is masked from 'package:purrr':
## 
##     simplify
```

```
## The following objects are masked from 'package:base':
## 
##     aperm, apply, rowsum
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'XVector'
```

```
## The following object is masked from 'package:purrr':
## 
##     compact
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:DelayedArray':
## 
##     type
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: bumphunter
```

```
## Loading required package: foreach
```

```
## 
## Attaching package: 'foreach'
```

```
## The following objects are masked from 'package:purrr':
## 
##     accumulate, when
```

```
## Loading required package: iterators
```

```
## Loading required package: locfit
```

```
## locfit 1.5-9.1 	 2013-03-22
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```r
mvals <- getM(mset_noob)

probe_anno <- readRDS('../../data/main/interim/1_1_probe_anno.rds')

# color key
pheatmap_color_code <- readRDS('../../data/main/interim/1_1_color_code.rds')

color_code <- readRDS('../../data/main/interim/2_3_color_code.rds')
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)

# DMCs
dmcs <- readRDS('../../data/main/interim/2_4_dmcs.rds')

# annotation
anno <- readRDS('Z:/Victor/Repositories/EPIC_annotation/hg19_epic_annotation.rds')

#color code
color_code <- readRDS('../../data/main/interim/2_3_color_code.rds')
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)

# 450k annotation
anno_450k <- read_csv('Z:/Victor/Data/DNAm annotations/HumanMethylation450_15017482_v1-2.csv',
                      skip = 7)
```

```
## Parsed with column specification:
## cols(
##   .default = col_character(),
##   AddressA_ID = col_double(),
##   AddressB_ID = col_double(),
##   Genome_Build = col_double(),
##   MAPINFO = col_double(),
##   Coordinate_36 = col_double(),
##   Random_Loci = col_logical(),
##   Methyl27_Loci = col_logical(),
##   Enhancer = col_logical(),
##   DHS = col_logical()
## )
```

```
## See spec(...) for full column specifications.
```

```
## Warning: 1739 parsing failures.
##    row           col expected actual                                                                    file
##   7101 Coordinate_36 a double  MULTI 'Z:/Victor/Data/DNAm annotations/HumanMethylation450_15017482_v1-2.csv'
##  10137 Coordinate_36 a double  MULTI 'Z:/Victor/Data/DNAm annotations/HumanMethylation450_15017482_v1-2.csv'
##  52780 Coordinate_36 a double  MULTI 'Z:/Victor/Data/DNAm annotations/HumanMethylation450_15017482_v1-2.csv'
## 120983 Coordinate_36 a double  MULTI 'Z:/Victor/Data/DNAm annotations/HumanMethylation450_15017482_v1-2.csv'
## 121831 Coordinate_36 a double  MULTI 'Z:/Victor/Data/DNAm annotations/HumanMethylation450_15017482_v1-2.csv'
## ...... ............. ........ ...... .......................................................................
## See problems(...) for more details.
```

## Remove samples

remove contamined and non-interesting samples


```r
pDat_filt <- pDat %>% 
  filter(maternal_contamination_norm_flip < 0.35,
         !Sample_Name %in% c('PM364_hofb_cs', 'PL293_v_R2', 'PM366_vc_R2', 'P131_hofb_cs'),
         !Tissue %in% c('Villi maternal', 'Trophoblasts enz', 'Mixture cs', 'Dead Cells and Lymphocytes cs'))

# filter to first trimester
mvals_filt <- mvals[rownames(betas),pDat_filt$Sample_Name]
betas_filt <- betas[,pDat_filt$Sample_Name]
```

# PMD

## Filter cpgs

D.I schroeder state that they remove probes in  the following regions:

* promoters
* cpg islands
* cpg island shores



```r
# find cpgs in pmds and meeting above criteria
anno_pmd_filt <- anno %>%
  filter(!cpg_id %in% c('island', 'shore'),
         !grepl('promoter', genes_id),
         !is.na(pmd_id))
nrow(anno); nrow(anno_pmd_filt) # 867052; 133370
```

```
## [1] 867052
```

```
## [1] 133370
```

```r
pmd_cpgs_epic <- intersect(anno_pmd_filt$cpg, rownames(betas_filt)) # 111378
pmd_cpgs_450k <- intersect(anno_450k$IlmnID, pmd_cpgs_epic)

length(pmd_cpgs_epic);length(pmd_cpgs_450k) # 111378;39142
```

```
## [1] 111378
```

```
## [1] 39142
```

## Calculate densities

Now I can calculate the density per celltype / trimeseter. I follow the same code from my enrichment
analysis.


```r
densities <- pDat_filt %>%
  
  # get sample names for each group
  select(Trimester, Tissue, Sample_Name) %>%
  dplyr::rename(Celltype = Tissue) %>%
  group_by(Trimester, Celltype) %>%
  summarize(Samples =list(Sample_Name)) %>%
  
  # duplicate rows
  mutate(array = 'EPIC') %>%
  ungroup() %>%
  bind_rows( identity(.) %>% mutate(array = '450k')) %>%
  
  # Calculate densities
  mutate(densities = case_when(
    array == 'EPIC' ~ map(Samples, ~as.vector(betas_filt[pmd_cpgs_epic, .x]) %>% density()),
    array == '450k' ~map(Samples, ~as.vector(betas_filt[pmd_cpgs_450k, .x]) %>% density())
    ),
         x = map(densities, 'x'),
         y = map(densities, 'y')) %>%
     
 # clean up results
 # remove input data
  select(Trimester, Celltype, contains('x'), contains('y')) %>%
  unnest()
```

## Plot

density plots


```r
### EPIC
densities %>% 
  filter(array == 'EPIC') %>%
  ggplot(data = .) +
  geom_line(size = 2,aes(x = x, y = y, color = Celltype)) + theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_grid(Trimester~Celltype) +
  scale_color_manual(values= color_code_tissue[unique(densities$Celltype)]) +
  labs(x = '% methylation', y = 'density', color = '', title = 'EPIC')
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
### 450k
densities %>% 
  filter(array == '450k') %>%
  ggplot(data = .) +
  geom_line(size = 2,aes(x = x, y = y, color = Celltype)) + theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_grid(Trimester~Celltype) +
  scale_color_manual(values= color_code_tissue[unique(densities$Celltype)]) +
  labs(x = '% methylation', y = 'density', color = '', title = '450k')
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

Let's look at a few specific PMDs

subset to those with a lot of cpg coverage


```r
# manually explore regions to plot
anno %>%
  
  filter(cpg %in% rownames(betas_filt)) %>%
  
  # count number of cpgs per pmd, parse the pmd_id into position
  group_by(pmd_id) %>%
  summarize(n = n(), pmd_width = paste(unique(pmd_width), collapse = ', ')) %>%
  filter(!is.na(pmd_id)) %>%
  separate(pmd_id, into = c('chr', 'start', 'end')) %>%
  mutate(chr = factor(chr, levels = paste0('chr', c(1:22, 'X'))),
         start = as.numeric(start),
         end = as.numeric(end)) %>%
  arrange(chr, start) %>%
  
  # filter out pmds with low coverage
  filter(n > 20) %>%
  datatable()
```

```
## Warning: Expected 3 pieces. Additional pieces discarded in 18 rows [396,
## 1285, 1448, 1564, 1714, 2177, 2292, 2542, 3411, 4079, 4104, 4110, 4378,
## 4654, 4709, 4867, 4870, 5162].
```

<!--html_preserve--><div id="htmlwidget-0c59e217fa553c5ef20b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-0c59e217fa553c5ef20b">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283","284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299","300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315","316","317","318","319","320","321","322","323","324","325","326","327","328","329","330","331","332","333","334","335","336","337","338","339","340","341","342","343","344","345","346","347","348","349","350","351","352","353","354","355","356","357","358","359","360","361","362","363","364","365","366","367","368","369","370","371","372","373","374","375","376","377","378","379","380","381","382","383","384","385","386","387","388","389","390","391","392","393","394","395","396","397","398","399","400","401","402","403","404","405","406","407","408","409","410","411","412","413","414","415","416","417","418","419","420","421","422","423","424","425","426","427","428","429","430","431","432","433","434","435","436","437","438","439","440","441","442","443","444","445","446","447","448","449","450","451","452","453","454","455","456","457","458","459","460","461","462","463","464","465","466","467","468","469","470","471","472","473","474","475","476","477","478","479","480","481","482","483","484","485","486","487","488","489","490","491","492","493","494","495","496","497","498","499","500","501","502","503","504","505","506","507","508","509","510","511","512","513","514","515","516","517","518","519","520","521","522","523","524","525","526","527","528","529","530","531","532","533","534","535","536","537","538","539","540","541","542","543","544","545","546","547","548","549","550","551","552","553","554","555","556","557","558","559","560","561","562","563","564","565","566","567","568","569","570","571","572","573","574","575","576","577","578","579","580","581","582","583","584","585","586","587","588","589","590","591","592","593","594","595","596","597","598","599","600","601","602","603","604","605","606","607","608","609","610","611","612","613","614","615","616","617","618","619","620","621","622","623","624","625","626","627","628","629","630","631","632","633","634","635","636","637","638","639","640","641","642","643","644","645","646","647","648","649","650","651","652","653","654","655","656","657","658","659","660","661","662","663","664","665","666","667","668","669","670","671","672","673","674","675","676","677","678","679","680","681","682","683","684","685","686","687","688","689","690","691","692","693","694","695","696","697","698","699","700","701","702","703","704","705","706","707","708","709","710","711","712","713","714","715","716","717","718","719","720","721","722","723","724","725","726","727","728","729","730","731","732","733","734","735","736","737","738","739","740","741","742","743","744","745","746","747","748","749","750","751","752","753","754","755","756","757","758","759","760","761","762","763","764","765","766","767","768","769","770","771","772","773","774","775","776","777","778","779","780","781","782","783","784","785","786","787","788","789","790","791","792","793","794","795","796","797","798","799","800","801","802","803","804","805","806","807","808","809","810","811","812","813","814","815","816","817","818","819","820","821","822","823","824","825","826","827","828","829","830","831","832","833","834","835","836","837","838","839","840","841","842","843","844","845","846","847","848","849","850","851","852","853","854","855","856","857","858","859","860","861","862","863","864","865","866","867","868","869","870","871","872","873","874","875","876","877","878","879","880","881","882","883","884","885","886","887","888","889","890","891","892","893","894","895","896","897","898","899","900","901","902","903","904","905","906","907","908","909","910","911","912","913","914","915","916","917","918","919","920","921","922","923","924","925","926","927","928","929","930","931","932","933","934","935","936","937","938","939","940","941","942","943","944","945","946","947","948","949","950","951","952","953","954","955","956","957","958","959","960","961","962","963","964","965","966","967","968","969","970","971","972","973","974","975","976","977","978","979","980","981","982","983","984","985","986","987","988","989","990","991","992","993","994","995","996","997","998","999","1000","1001","1002","1003","1004","1005","1006","1007","1008","1009","1010","1011","1012","1013","1014","1015","1016","1017","1018","1019","1020","1021","1022","1023","1024","1025","1026","1027","1028","1029","1030","1031","1032","1033","1034","1035","1036","1037","1038","1039","1040","1041","1042","1043","1044","1045","1046","1047","1048","1049","1050","1051","1052","1053","1054","1055","1056","1057","1058","1059","1060","1061","1062","1063","1064","1065","1066","1067","1068","1069","1070","1071","1072","1073","1074","1075","1076","1077","1078","1079","1080","1081","1082","1083","1084","1085","1086","1087","1088","1089","1090","1091","1092","1093","1094","1095","1096","1097","1098","1099","1100","1101","1102","1103","1104","1105","1106","1107","1108","1109","1110","1111","1112","1113","1114","1115","1116","1117","1118","1119","1120","1121","1122","1123","1124","1125","1126","1127","1128","1129","1130","1131","1132","1133","1134","1135","1136","1137","1138","1139","1140","1141","1142","1143","1144","1145","1146","1147","1148","1149","1150","1151","1152","1153","1154","1155","1156","1157","1158","1159","1160","1161","1162","1163","1164","1165","1166","1167","1168","1169","1170","1171","1172","1173","1174","1175","1176","1177","1178","1179","1180","1181","1182","1183","1184","1185","1186","1187","1188","1189","1190","1191","1192","1193","1194","1195","1196","1197","1198","1199","1200","1201","1202","1203","1204","1205","1206","1207","1208","1209","1210","1211","1212","1213","1214","1215","1216","1217","1218","1219","1220","1221","1222","1223","1224","1225","1226","1227","1228","1229","1230","1231","1232","1233","1234","1235","1236","1237","1238","1239","1240","1241","1242","1243","1244","1245","1246","1247","1248","1249","1250","1251","1252","1253","1254","1255","1256","1257","1258","1259","1260","1261","1262","1263","1264","1265","1266","1267","1268","1269","1270","1271","1272","1273","1274","1275","1276","1277","1278","1279","1280","1281","1282","1283","1284","1285","1286","1287","1288","1289","1290","1291","1292","1293","1294","1295","1296","1297","1298","1299","1300","1301","1302","1303","1304","1305","1306","1307","1308","1309","1310","1311","1312","1313","1314","1315","1316","1317","1318","1319","1320","1321","1322","1323","1324","1325","1326","1327","1328","1329","1330","1331","1332","1333","1334","1335","1336","1337","1338","1339","1340","1341","1342","1343","1344","1345","1346","1347","1348","1349","1350","1351","1352","1353","1354","1355","1356","1357","1358","1359","1360","1361","1362","1363","1364","1365","1366","1367","1368","1369","1370","1371","1372","1373","1374","1375","1376","1377","1378","1379","1380","1381","1382","1383","1384","1385","1386","1387","1388","1389","1390","1391","1392","1393","1394","1395","1396","1397","1398","1399","1400","1401","1402","1403","1404","1405","1406","1407","1408","1409","1410","1411","1412","1413","1414","1415","1416","1417","1418","1419","1420","1421","1422","1423","1424","1425","1426","1427","1428","1429","1430","1431","1432","1433","1434","1435","1436","1437","1438","1439","1440","1441","1442","1443","1444","1445","1446","1447","1448","1449","1450","1451","1452","1453","1454","1455","1456","1457","1458","1459","1460","1461","1462","1463","1464","1465","1466","1467","1468","1469","1470","1471","1472","1473","1474","1475","1476","1477","1478","1479","1480","1481","1482","1483","1484","1485","1486","1487","1488","1489","1490","1491","1492","1493","1494","1495","1496","1497","1498","1499","1500","1501","1502","1503","1504","1505","1506","1507","1508","1509","1510","1511","1512","1513","1514","1515","1516","1517","1518","1519","1520","1521","1522","1523","1524","1525","1526","1527","1528","1529","1530","1531","1532","1533","1534","1535","1536","1537","1538","1539","1540","1541","1542","1543","1544","1545","1546","1547","1548","1549","1550","1551","1552","1553","1554","1555","1556","1557","1558","1559","1560","1561","1562","1563","1564","1565","1566","1567","1568","1569","1570","1571","1572","1573","1574","1575","1576","1577","1578","1579","1580","1581","1582","1583","1584","1585","1586","1587","1588","1589","1590","1591","1592","1593","1594","1595","1596","1597","1598","1599","1600","1601","1602","1603","1604","1605","1606","1607","1608","1609","1610","1611","1612","1613","1614","1615","1616","1617","1618","1619","1620","1621","1622","1623","1624","1625","1626","1627","1628","1629","1630","1631","1632","1633","1634","1635","1636","1637","1638","1639","1640","1641","1642","1643","1644","1645","1646","1647","1648","1649","1650","1651","1652","1653","1654","1655","1656","1657","1658","1659","1660","1661","1662","1663","1664","1665","1666","1667","1668","1669","1670","1671","1672","1673","1674","1675","1676","1677","1678","1679","1680","1681","1682","1683","1684","1685","1686","1687","1688","1689","1690","1691","1692","1693","1694","1695","1696","1697","1698","1699","1700","1701","1702","1703","1704","1705","1706","1707","1708","1709","1710","1711","1712","1713","1714","1715","1716","1717","1718","1719","1720","1721","1722","1723","1724","1725","1726","1727","1728","1729","1730","1731","1732","1733","1734","1735","1736","1737","1738","1739","1740","1741","1742","1743","1744","1745","1746","1747","1748","1749","1750","1751","1752","1753","1754","1755","1756","1757","1758","1759","1760","1761","1762","1763","1764","1765","1766","1767","1768","1769","1770","1771","1772","1773","1774","1775","1776","1777","1778","1779","1780","1781","1782","1783","1784","1785","1786","1787","1788","1789","1790","1791","1792","1793","1794","1795","1796","1797","1798","1799","1800","1801","1802","1803","1804","1805","1806","1807","1808","1809","1810","1811","1812","1813","1814","1815","1816","1817","1818","1819","1820","1821","1822","1823","1824","1825","1826","1827","1828","1829","1830","1831","1832","1833","1834","1835","1836","1837","1838","1839","1840","1841","1842","1843","1844","1845","1846","1847","1848","1849","1850","1851","1852","1853","1854","1855","1856","1857","1858","1859","1860","1861","1862","1863","1864","1865","1866","1867","1868","1869","1870","1871","1872","1873","1874","1875","1876","1877","1878","1879","1880","1881","1882","1883","1884","1885","1886","1887","1888","1889","1890","1891","1892","1893","1894","1895","1896","1897","1898","1899","1900","1901","1902","1903","1904","1905","1906","1907","1908","1909","1910","1911","1912","1913","1914","1915","1916","1917","1918","1919","1920","1921","1922","1923","1924","1925","1926","1927","1928","1929","1930","1931","1932","1933","1934","1935","1936","1937","1938","1939","1940","1941","1942","1943","1944","1945","1946","1947","1948","1949","1950","1951","1952","1953","1954","1955","1956","1957","1958","1959","1960","1961","1962","1963","1964","1965","1966","1967","1968","1969","1970","1971","1972","1973","1974","1975","1976","1977","1978","1979","1980","1981","1982","1983","1984","1985","1986","1987","1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","2023","2024","2025","2026","2027","2028","2029","2030","2031","2032","2033","2034","2035","2036","2037","2038","2039","2040","2041","2042","2043","2044","2045","2046","2047","2048","2049","2050","2051","2052","2053","2054","2055","2056","2057","2058","2059","2060","2061","2062","2063","2064","2065","2066","2067","2068","2069","2070","2071","2072","2073","2074","2075","2076","2077","2078","2079","2080","2081","2082","2083","2084","2085","2086","2087","2088","2089","2090","2091","2092","2093","2094","2095","2096","2097","2098","2099","2100","2101","2102","2103","2104","2105","2106","2107","2108","2109","2110","2111","2112","2113","2114","2115","2116","2117","2118","2119","2120","2121","2122","2123","2124","2125","2126","2127","2128","2129","2130","2131","2132","2133","2134","2135","2136","2137","2138","2139","2140","2141","2142","2143","2144","2145","2146","2147","2148","2149","2150","2151","2152","2153","2154","2155","2156","2157","2158","2159","2160","2161","2162","2163","2164","2165","2166","2167","2168","2169","2170","2171","2172","2173","2174","2175","2176","2177","2178","2179","2180","2181","2182","2183","2184","2185","2186","2187","2188","2189","2190","2191","2192","2193","2194","2195","2196","2197","2198","2199","2200","2201","2202","2203","2204","2205","2206","2207","2208","2209","2210","2211","2212","2213","2214","2215","2216","2217","2218","2219","2220","2221","2222","2223","2224","2225","2226","2227","2228","2229","2230","2231","2232","2233","2234","2235","2236","2237","2238","2239","2240","2241","2242","2243","2244","2245","2246","2247","2248","2249","2250","2251","2252","2253","2254","2255"],["chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr2","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr3","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr4","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr5","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr6","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr7","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr8","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr9","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr10","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr11","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr12","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr13","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr15","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr16","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr17","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr18","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr19","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr20","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr21","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22","chr22"],[1042309,1062711,1090448,1100009,1105870,2336017,2405367,2468224,2517364,2707094,2767222,2818155,2840734,2876784,2928904,2980579,2989839,3028204,3049129,3071153,3101770,3154504,3194314,3230122,3265145,3300896,3341493,3597414,3613382,3895129,4121533,4752106,5759967,6034424,6054573,6131627,6944896,7045315,7688444,11325385,11519210,11842143,12748427,14033417,14092325,18000001,18025430,18310261,18576084,18682255,18845892,18916482,18983520,20500001,20566122,20707471,29530869,29800670,30670294,33750001,33871707,34415612,36779762,37273212,38003476,38321013,38817121,47110360,47688540,47963580,48235564,48710271,49137458,50286908,50572124,57054617,57663226,60250001,60312129,68735595,69807084,70750001,71319372,72522325,74436624,75006722,76062367,76313242,77107123,78500001,78730105,79245126,90288363,95505047,97255529,98292572,99243557,100778474,101617277,103893839,106425378,107485987,108617396,110867761,110899886,110951710,111019461,111124957,111308531,112250001,112335297,114405420,114498709,115199466,144926008,146614717,147034073,147271085,147781217,150286389,150455698,150754895,150999963,151643632,154913885,155431760,155820249,156215678,156386329,156418128,157425123,159767377,159849002,160220664,160306677,161250001,161659725,163170133,164403073,166863195,168310622,168900803,173693002,173835432,174500001,175417789,179500001,179554497,179649472,179721885,180034517,184750001,184916705,186656886,191664775,201635769,202921509,203064554,203217834,206828364,208178678,209091370,213000001,213110627,213323260,214018414,215102278,216165289,217007923,217853186,219135783,220298763,228000001,228500001,228559569,228629326,231603460,231816938,232108081,236041802,236526898,239053533,240755306,245823463,245869765,246167350,246922094,36560,723097,855092,937287,1042731,1153821,1265650,1473247,1495381,1575028,1649211,1822145,1874940,1920568,2172724,2633366,2752626,3042680,3677341,3729803,3812194,4134320,5083317,6553612,7489094,7825192,8085738,10506607,10584217,10971404,11018940,13075051,17924402,19426447,21120467,29314632,30792700,31660287,34172243,39989527,40632004,44883540,45095084,48950440,50428322,51109003,52710629,59382383,60666259,67814056,75837591,76964300,79074101,79593752,80403310,81548193,85808173,87011521,87023815,89082415,95556800,95678958,95819573,97296848,97573657,103000001,106575293,106870649,107969900,109113266,109828018,109876203,109909995,113912034,113978265,114028661,114500001,115136732,116353725,118861532,119333297,119485867,121061812,121128786,121463322,122252384,124499726,127360575,129211545,129727304,131071702,131774567,131838233,132891915,133145341,134083414,135000001,135193464,136593522,136897772,137240167,139255121,139376866,142605024,143462892,144998740,148030898,150018014,151250001,153430988,154043705,154438165,155264316,157250001,158069457,158594022,162020918,166554899,166940958,175771001,177263950,178192030,180580120,181553760,182031275,184250001,185171835,192769058,194387198,198291033,198359845,198378400,207017153,207215668,209250001,211250001,212468910,214022502,217267730,217383460,219482348,219576299,219905013,220086191,220923195,223521995,223626282,225615709,226155976,228500001,228754798,234147538,234328811,234512597,236811316,237472409,239105466,239207293,240000001,240114702,240231582,240315143,240403171,240630539,241261351,241290033,241370788,241507209,241554585,242601511,242638531,242680236,215141,325894,2116785,3196618,6879642,10345807,10375587,10522390,10833448,11034052,11863539,12021628,13221040,14572656,18462061,19165105,21094267,22389370,26641169,27746947,28365860,35656633,36397215,38666518,39519207,39826825,42282525,43774183,51725102,62339354,62836228,63239243,64774715,67006402,67131961,68000001,68139855,70099353,73667899,73757063,76038996,77230087,79712800,82939932,85091525,87221234,87359530,87925251,88345718,89246795,96139847,98016313,99297444,100250001,102000001,103500001,114044436,115354501,115553346,116349812,116986269,119199252,122000001,122110267,125235887,128337536,128539984,128657632,132270166,132728709,133236827,135998668,138405524,140246801,140886646,141137803,145174967,147854493,148625079,149898620,154750001,158743916,160036817,162632622,165204987,167271816,168580890,168935522,169750001,170350105,170869879,174076249,177186162,181237937,182190319,183156569,187562803,188400448,192063591,192646812,193250001,193610686,194265725,194442065,609913,1039924,1097729,1156583,1399225,1552172,1571021,3511174,3536984,3655338,3763069,4831227,4925075,5104994,5946018,6061014,6165602,6253626,6431447,6525436,7246369,7563329,7644487,7750001,8052688,8149975,8211794,10213018,11040355,13238754,16250001,18297687,19865967,20426079,21559964,22299966,26639968,30333941,31492167,32527753,35922785,42095560,42354308,44423477,46087330,46729698,47403470,48750001,48858514,52638006,57918320,59479327,62065844,63898180,66219270,68500001,68898517,69375442,70950653,72839998,73116829,73654186,74954409,76204129,78421692,81526977,82355560,85086705,86250001,92258714,93446031,94975334,96754986,101355130,104860247,107527873,109313996,111780258,114216630,115120260,115740382,117846559,120887994,121207605,127506395,130293405,133118641,134627906,138045649,147796213,148086610,150332976,154933691,155474054,155632844,155691133,155922779,156517545,160666716,164307597,164750001,165524628,168008608,171442051,172971694,174754559,176692757,177550450,177951452,178849744,181732835,188771908,189569997,189817757,190522529,190804717,63001,590523,618968,838710,980402,1012173,1063040,1078630,1091427,1165389,1209467,1223942,1398470,1439863,1470005,1499649,1721082,1800051,2091950,2165485,2190306,2233468,2617753,2696310,2920006,3241376,3379646,3589356,3660136,4015513,4919631,5940529,6431847,6502630,8000001,8511274,9391035,9599716,11000001,11438522,13000001,13863334,15150038,15553924,15990317,16233421,17688473,20341694,21838219,24309719,25226871,25945409,26705662,28964395,32750208,33589611,33928161,33974067,37827825,37876484,41360719,42000001,43803431,44845450,50301495,50731110,52131783,52321247,52441525,54250001,57868241,59225652,61000001,62507802,63497863,63838362,70399498,71003580,71051472,72714178,75098966,76824199,81750001,82750001,82806356,84000001,86804071,87934493,88021567,88221346,89890713,90800422,91711885,94114895,95567484,95794837,96750001,97612887,100561449,110127513,110456170,110588346,113000001,113419905,113727095,116729999,118981758,119829094,121545881,125028936,129118170,129269512,132899699,134822192,134855544,134907719,134942881,135216410,135294417,135392936,135444375,135670223,136144755,140169171,140244339,140327116,140462645,140561810,140607558,143565258,143958678,145195333,145296548,145706046,145739076,146329292,146594754,146814132,146870034,147545030,148014274,151285018,154390914,155123021,157363411,158466579,158623120,158691418,159824527,160907963,163380517,166338496,167218446,168660854,168997281,169686981,169864050,173349193,174346042,175018363,175157286,175232673,176007230,176040128,176103403,180250001,180475009,338791,771324,832909,1076448,1262963,5954299,6003664,17077646,18639748,22678685,27969575,28890282,29241953,29629813,29708641,29725954,29940066,29964545,30003097,30019347,30032513,30053149,30083349,30137175,30151480,30179408,30203590,30500001,30527610,31010426,31031765,31273903,31348100,31384506,31433191,31479897,31574168,32272362,32660310,32741006,32837856,32930349,32970764,33083905,33156793,33192974,33237697,33269316,33280802,33284366,33289235,34000001,39024124,40500001,40675559,41267474,41450522,46567217,48701895,49750001,50800335,51529382,54083836,55147352,62236549,64899080,66068203,66861553,69402205,70750001,71049634,72653126,72949599,73389831,74569442,75851809,75973108,76759896,78417064,83132039,84278060,84476098,84994448,85541438,86445878,87704427,88814599,88934169,91377941,94186399,95937265,96132764,96570858,97076822,97837921,99502577,100173999,100548874,100783820,101523653,101957696,103476330,115025740,117500001,118348263,123250001,123359453,124167009,125726565,126703081,127006682,127483454,130729430,133315037,133605599,134027154,134201043,134258974,146327395,146392474,146797595,146961939,152123913,152807712,153000166,153524294,155677575,156760240,159511308,159575388,159656964,160500224,160599692,160704126,161832952,162334772,164750001,166188452,166342232,166503414,166789221,166961670,167507550,167554131,167685115,169500001,170763701,60976,126515,177676,215842,915695,961101,1254185,1282474,1337390,1390580,1576564,1627308,1677044,1720233,1775887,2500001,3056533,4823693,4842756,4868280,4890103,5089691,8450237,9732741,14152359,18585721,19151559,19779531,25250001,27137164,27180922,27186276,27195569,27199623,27206837,27211823,31343068,32077436,32305251,36500001,38467957,38637534,45164333,46364354,48465300,48934895,49786299,50169533,51500001,53677655,54582292,63290398,66500001,70000001,70236319,75669479,80750001,81078645,81911457,82630349,84407504,84654179,86254049,86527175,87886089,88227220,90064301,92512276,92830464,93042318,95240369,97938226,98085696,100396514,101922581,102044710,102099569,102874077,103757196,108076210,112367454,117371618,117641958,119702722,122346651,124357436,126820321,131000001,131281037,131913068,135436293,140361413,140421030,141235858,141969929,143000001,143028127,143524043,144258409,145445018,149250001,149377403,149776496,150048952,150129140,150185675,152254239,153381378,154001568,154494179,155316656,155367703,155549111,155591148,155608475,155641416,155957783,157062156,157080084,157196906,157219503,157312064,157324377,157340769,157362857,157396047,157650035,157689789,157754508,157803643,157913471,157943259,158750001,780799,840406,979270,1041212,1959814,1991467,2065835,2468328,3026776,4840044,4968355,9750001,9828959,10376885,10735045,10911585,10966622,11096653,11243100,11362108,11508659,11604622,12567598,13286346,15675114,16687928,19792356,19842287,20295614,23542234,23640707,25589682,26578013,26780484,31161241,31617889,32526444,33700143,35213654,40521359,47648540,49250001,49456101,50053428,50985414,52484775,53250001,54016980,55696587,56178316,56524723,57331654,57395350,57523370,58069007,58355439,62214986,62951259,63324619,64500001,65453501,65662700,65874277,68595585,69027501,70124629,70723015,72430344,72632116,72919222,73326746,73612280,75195599,75425600,76059672,76482752,77757416,78075613,79879959,80688321,82917032,85259803,86909966,87750001,89000001,89409598,92375253,96928560,97243057,98359581,104583090,105305156,105449163,105670881,106401297,108987317,109165026,109869024,110773074,111057249,114514262,116750798,118155409,118602493,119193614,119703495,120033819,123291796,125044319,126780930,127816588,131775958,132123932,133133139,133210999,133562769,133757181,134927216,136250001,136772282,138245924,139578957,139959630,140700515,142521693,142598212,142840065,142982561,143052956,143200776,143543181,143693169,143913408,144106277,144313489,144339200,7266158,8848709,11608589,12765814,17569805,17897489,18428803,22437696,23812668,23841523,25750001,27958560,29816737,31244764,32563228,37028342,71642206,71849098,72225269,74983835,77219085,78500001,80259024,82166233,83983054,86760959,95500001,95663303,99606001,102500001,102831995,103540898,107750001,115517434,116892869,117956997,118760419,119016832,119217318,119547464,121172036,123961992,124250001,127817800,132068385,133945157,134275719,134352303,135441098,135464570,136252230,136498624,136533127,136674219,136800564,136874297,137208116,137296805,137333627,137572547,139897229,139936610,140037874,140061348,140134692,140164935,140190013,140231412,1221215,1242506,1274423,1406250,1430015,1507314,1585793,1770019,2347472,2555223,3000001,5113574,6443766,6662649,9490216,10348596,11100531,11247987,13974191,17000001,18064682,18470241,19356165,20342154,23528550,23807247,24024985,26152267,26264029,31500001,43750001,43821277,44250001,45641597,45896971,46061237,48277916,49483359,49534608,50010185,50492673,50557869,50857417,50895224,51498294,52170028,54500001,57061222,57791282,60792543,63750001,63804427,64248785,65250001,65403583,65600897,66484068,71232883,76997116,77542312,79510988,79733595,81000001,81130750,82500001,83625710,85314983,85889621,85945284,85987408,86500001,87829696,88013260,90000001,90333202,90630604,96442781,98856371,99600107,99781311,100289257,100984179,101080646,102312750,102409659,102499637,102580393,102890254,102896684,102986637,106088547,106392803,108914796,109664955,110216456,110662317,112829294,116844085,117875020,118024221,118890320,118925171,119303554,119484982,119912422,120345812,121816274,122207368,124392955,124568434,124629784,125000001,127000001,129233298,129427357,129899611,129954269,130229768,130398649,131482607,131556306,131587403,131661687,132000001,132255972,132474155,132870687,133032833,133208315,133308486,133358582,133494805,133729246,133788148,133869767,133891775,133922184,133945811,134474266,134532908,134548959,134591024,134606843,134657436,134794783,134805791,134831050,134852974,135192552,135240969,1026828,1046742,1073359,1086520,1098042,1288075,1316983,1388411,1552720,1630801,1672579,1750001,2013893,2144822,2267083,2393228,2398979,2423395,3098454,3143345,3181654,3211053,4250001,4586107,5915703,6633673,6661520,6904789,6998526,7230812,8366955,13768595,13942926,14952485,15093122,17455203,17522629,17714863,18084288,18770656,20142074,20579976,20648422,22171897,22319954,22411095,22645194,22807944,23708938,24475315,26310389,26677593,27701141,28504885,29996871,30301523,30564705,31805353,33491242,36262221,36488492,43559792,44298334,44706338,45000001,45264552,45333820,45644072,48156691,55398697,57597345,58250001,58808174,59422434,60039622,62547922,79059339,83221432,86750001,88000001,89160758,89947092,95783522,95942795,98397032,100505648,102967414,103540269,104986633,105877463,108523113,108800194,110113614,111750001,112339701,113270432,113575328,113910835,114337265,115631760,116894682,117173352,117529405,117735579,118757775,118799154,118961313,120023971,120362385,120400293,120534138,120669163,120739206,120828513,123030896,123117756,124541637,125064593,125914276,126379020,130437625,131066824,131286743,132458518,132736856,133208919,133332668,133762842,133788130,133841713,134000001,134087792,134111074,134138198,134337374,1888009,4750001,4891433,7269718,9250001,11250001,13316756,13608302,14026510,18626881,19750001,19985626,20414390,20597272,21753044,22931113,24607511,28632309,29828701,33484318,34394882,38905181,39373370,39869243,40253743,43557103,46098993,47042254,47208397,50271552,51281580,51394739,52431738,53268859,56313324,56355042,56750001,57790701,58276775,66500001,69290325,70626753,70953819,76000001,76243572,79305929,79855646,79996251,80676806,82126755,84767044,95500001,95825419,96545248,97813506,98902714,99635754,101517784,101742819,101883760,102220549,102413974,106059007,106693701,106821874,107047983,109914759,109959626,111497913,112402491,112514136,112888761,113370963,116283832,116683671,117544828,117696777,118256978,124238066,124381771,124584319,124734612,124912485,125242149,125778483,126331469,126506861,126811790,127466019,127995797,128354114,128750806,129331501,129390037,129519783,130215562,130350153,130484039,130509715,130543797,130611338,131501061,26030678,27269769,28292899,28812302,29068247,30146619,30379106,30673035,33739854,34383359,35250001,35904454,37164169,37342229,39250001,46349434,52212047,52323997,52674347,54557080,57106932,60056450,63316681,66703461,68457306,69580221,71338061,73607603,75360322,77543259,78131419,79916831,80127447,83351899,85708326,87127887,88813777,90849718,92678879,95542164,100584959,101250001,101367497,102500001,106500001,107319065,107669330,107947021,110862024,111036895,111073796,111119131,111551996,111809115,111897223,111943710,112046276,112308927,112600775,112719193,112753780,113616115,113766854,113804640,113859683,113945156,113983943,21156549,22929461,22959249,23981960,24589453,26136419,28575225,29467172,35517399,36206100,36737245,37161928,37795289,39000001,41146080,44859130,46190910,47215340,51606473,51851720,55432690,55655988,56354312,57000001,57402875,59069604,59167307,59264056,59864621,59934431,60052390,61654033,65495745,69582748,69725831,75750001,75889133,77517478,79500001,81070188,85070232,91750001,91860466,92967220,93325487,93711402,95250602,95413324,96500001,96569930,96755700,97628590,100104983,100123178,100191119,100271517,100364187,100419614,100602138,100995749,101123917,101171350,103563125,103760028,103872023,103968655,104011512,104142509,105259906,105306839,105404966,105741746,105814883,105928833,106061378,106292344,22752148,24000001,24799603,25270607,25460656,25922762,26322106,28183871,28305738,29750001,30110353,30247957,31347127,31391296,35251301,35570414,43784169,50885769,51839297,52238401,57944364,59146849,59309009,67679088,70876750,71132423,71398583,71449069,71522663,75714418,75899909,77250001,77363326,77512699,78331718,78484449,79213758,80435233,80685622,80742720,80821720,81006719,82114776,84250001,86602022,89500001,89641345,90198687,90260683,90739526,95867175,95997396,96367673,98500001,100122673,696698,750591,815591,961305,1028782,1070396,1095321,1119944,1195083,1208607,1250908,2789147,2808004,3149415,3183585,5230388,6010403,6473985,7294610,7644092,8106570,9247181,12905192,14705157,17748357,19878009,19992807,20268000,21220835,21597754,21676127,22356530,22733961,22867378,23755604,25612133,27700111,28287476,31456686,31488525,34153611,34299876,34667376,47402766,49432824,49747762,51139112,52883205,53648449,53923524,54071028,54100776,54291960,54352441,54424500,57419988,58346574,60628136,69392800,69822161,69950251,70750001,72210877,74274802,74826695,77855722,78669224,80761733,81219315,84565189,84895396,84968959,85171323,85272827,85427647,85470172,86206190,86750001,86870055,86987241,87561795,87598406,87773462,88700427,88765725,31869,85415,171341,2899124,6000001,6299680,6501528,6558196,6620435,6738449,7970816,8750001,10123617,10653082,11086142,11443054,12000001,14657407,21161006,21260635,28404942,28644579,29391520,29508394,29989547,30594641,30811959,30839061,31535095,32240097,33368567,33555734,41124435,41218692,41330837,41750001,41806356,42023850,43166388,43211022,43222744,46794623,47591466,49189900,51278139,66281275,69962213,73099453,73301149,73395924,73467548,74645771,74705865,75224894,899575,3870088,4445648,5287100,5886250,6274616,6404874,10750001,11139937,11680109,12083843,13209045,13253272,13377687,13538810,13816524,13859027,20373259,21184903,22697457,23454272,24012174,28606301,29274846,29993521,30058874,30328291,30875688,32132486,33401627,35675828,39056543,41000001,42275679,43044677,47250001,48122758,51816324,52940193,55092521,55879886,59828768,61569621,63334995,64584769,66290530,71296909,71757647,72191836,72243780,72337627,72899579,72974415,73092811,73987790,74100112,74174547,74499367,74658298,74727758,74740733,74795184,1465035,1477682,2601912,3179993,3244549,3789224,4536173,4674982,5201669,5245212,5292062,6486080,6542073,7653715,8598199,8669795,9938445,12929153,12997170,13208347,13227159,13258646,13271132,13478268,14181645,14220973,14237962,14598167,14750001,15480614,15523359,15556536,15883936,17589648,17620411,17846134,24147611,32423623,33092674,33976851,34276105,34396172,34713208,35250001,35411811,35558331,35628467,35634022,35902725,36540050,38704777,38806194,38867609,39089920,39917157,39956116,40022022,40492872,40535028,40716644,40728448,40741514,41083430,45727281,46078430,46223892,46287981,46720542,47137096,47155536,47201935,47329525,47750001,47894154,49696687,51136049,51148344,51274285,54617808,55761413,55834799,55863836,55881922,55891160,56031323,56108062,56303989,56617226,1155120,2622699,2729498,4000001,4101399,4150766,4178497,5245799,6141876,6698958,8998019,9444894,9767862,9963915,10148385,10913072,12088872,16659116,17157623,18713701,20005790,21326246,21451932,21750001,22515881,23917434,24400038,30048440,30104448,31000001,34380500,35468618,35660256,36793268,36898262,37112516,39755374,40061440,41811568,41978552,42178438,42389143,42812853,42872961,43250001,43532951,44094356,45848915,46013575,46263486,46447378,47369092,50538708,52341119,52526401,52741955,55533657,55757661,56523644,57504117,57586243,57614183,58373506,59553108,59806935,60485293,60562349,60785943,60811666,60874039,61146263,61259115,61424854,61444742,61556389,62189461,62208516,16938187,18539746,22733938,27881303,29945844,30233978,33191225,34909669,34965095,35185558,35321130,36181473,38635646,44583954,44801958,45547236,45572745,45650663,45835592,45888737,45991355,46119645,46282168,46912535,19868850,19993798,20070182,20102443,20667334,20802745,20907913,21231699,21317215,23412113,23490213,23532554,23905543,24009274,24500001,24895987,25483151,26164631,30797169,30902526,31784506,32648638,35795278,35933700,35985786,37520066,37599107,37653475,38284615,38412337,41953612,42070230,42161107,42227734,42539746,42651847,42884839,43107959,43272714,43395153,43785037,45987987,46610314,46984172,47109057,47475475,47533215,47762665,48023894,48871926],[1062231,1083073,1099176,1105436,1122689,2349785,2423507,2507256,2527530,2761985,2817644,2836915,2866404,2908224,2969134,2987131,3018038,3046400,3061758,3092399,3137704,3172680,3229774,3264888,3296928,3311128,3360764,3613177,3624492,4092492,4613848,5335220,5840612,6054221,6110329,6161823,7045105,7645986,7716679,11461255,11631353,11893168,12975585,14091936,14314769,18024964,18306691,18564460,18681459,18829481,18915976,18983306,19053531,20542152,20683048,20751480,29677548,30600523,30930595,33863191,34222257,34981447,37270963,37556197,38029466,38714505,38946965,47261812,47746295,48222456,48593977,49014957,50262003,50571253,50653502,57660550,58487859,60311902,61062162,69805554,70157630,71284811,72521058,74436413,74938956,75363404,76250000,77105698,77334825,78729211,79244868,80559960,90944689,96174732,98158776,98894522,99502569,101133175,103636494,106424915,107400727,107780645,108765271,110899588,110950505,111017766,111124745,111307403,111422276,112325990,112656667,114496658,114687871,115370284,144988911,146727982,147120619,147404165,147842571,150346620,150754601,150999581,151532173,151714789,154930051,155735365,156215343,156230163,156414056,157239941,158000694,159842157,159908807,160259645,160580206,161376804,162449756,163353379,165074964,167341389,168897079,169370337,173834998,174034334,175400014,176136206,179553922,179648912,179717933,180034115,180250000,184915932,185309520,190710851,194844330,201690346,203064232,203217604,203292873,207460582,208384409,209373289,213110399,213322673,213759613,214462834,215329620,216404402,217399561,218141482,219787742,220705383,228264376,228535126,228627725,228844699,231815995,232106920,232416506,236329603,237616730,239586725,241120219,245869360,246068549,246921652,246974833,179812,853929,935685,1042314,1153015,1249017,1459764,1495178,1523278,1631001,1725839,1874703,1920074,2168206,2633114,2733748,3042378,3130513,3728702,3811689,4028523,4983316,5730932,6838857,7824780,8085472,8208029,10583991,10609139,11018481,11131775,14330522,18242553,19705790,23409823,29996334,30979818,31873304,36435900,40188589,41728572,45008698,45249372,50427548,51108025,52651906,53504120,60329295,60750000,68109544,76639223,79073699,79593203,80383187,81274714,84371224,85834009,87022683,87156032,89411551,95670570,95818899,96021212,97572551,97628793,104642019,106825954,107527132,108250133,109250000,109875359,109909086,110062340,113973244,114015504,114057770,115136197,115548520,118288581,119248534,119485612,119630595,121128452,121209326,121551503,124498721,127130165,127446141,129375768,130061512,131166642,131827318,131868618,133059583,133715888,134453047,135192167,135316264,136897530,137238929,138294954,139376655,142604193,143147086,144212392,147061270,148274309,150890468,151757319,154042694,154436151,155262250,156566357,157941497,158162044,159021417,162517662,166940633,167857278,176462095,177473081,178645840,181553035,182030005,182229465,185171452,186311636,193996935,195907746,198357995,198377323,199934015,207215018,207338122,209718337,211766711,213110424,215113362,217383209,217623940,219532394,219614350,219930771,220095318,222144277,223625819,224187296,226154959,226888509,228754200,229583652,234316984,234390187,234662670,236895764,237572655,239206976,239419844,240114010,240231110,240314934,240402963,240449069,240712373,241289081,241365711,241406813,241554376,241571449,242636313,242678972,242751149,325659,2115028,3040805,3511333,8425481,10375380,10522177,10832686,11009445,11170900,12020651,12153607,13298437,14590369,19163249,19454725,22388495,22876568,27329527,28245035,28591666,36396705,36720198,38948784,39826028,40061993,42437888,44011410,51936602,62834653,63239025,63504492,65282382,67131141,67498889,68139585,69063469,70984129,73755916,74746192,77171184,78337978,79898324,85090970,87122296,87358558,87924482,88146951,89246321,90587544,97978303,98959175,99706391,100958998,102523387,106555215,114175039,115552411,116348122,116859953,119198786,120089380,122109566,122617476,125470412,128539522,128657061,128730039,132563264,133236566,133428500,137064976,138962681,140484439,141136116,141441107,147242651,148559801,149898128,1,155278557,159294743,160543531,165204740,167120455,168580617,168935066,169204527,170346720,170862127,170965028,174595960,178006845,181802513,182895704,183430358,187750000,188870604,192407587,193042073,193608511,193714762,194441399,194635220,631534,1053740,1127997,1178184,1413146,1565251,1595074,3535108,3613566,3698791,3842889,4905513,5071903,5580955,6060768,6158030,6251743,6274466,6474113,6615651,7534441,7644268,7688588,7790167,8138068,8204976,8239378,10979326,12936984,14473605,17039596,19862373,20310914,21559182,21925454,23374883,30327950,31421166,32468752,35922025,36921644,42353034,44144142,46086659,46689884,47059018,47533360,48858213,49033698,52905416,58129209,61748477,63897441,66217787,67975021,68730656,69275441,70594206,71441665,73116427,73652646,73986438,75028732,76597007,78745572,82170964,82515049,85621787,86524109,93322639,94968517,95260909,99283624,101637594,105632103,108176003,109728338,113085104,114901517,115738870,117739878,118066534,121207064,121841014,128649534,132863796,134288611,136084829,138476798,148086278,148621647,151218699,155473579,155557345,155690582,155767433,156348617,156679728,164210796,164472403,165328832,165838073,168391579,172970308,173945577,175371833,177159466,177950236,178350278,181216031,183093876,188996119,189789378,190465553,190803542,190968639,103804,614825,655123,896457,1011784,1056936,1078329,1090963,1106543,1208222,1223701,1246722,1438977,1469791,1497258,1510170,1755284,1822313,2164835,2189796,2233250,2617496,2696104,2714314,3241098,3358586,3588629,3643643,3907535,4919316,5183285,6431290,6501752,6635881,8510431,9201589,9597691,10165084,11437680,11956549,13863125,14003445,15553075,15981087,16232063,16489487,18342531,21428297,23987551,25226259,25945098,26705304,28845156,31229708,33184864,33927947,33971924,34021279,37870427,38293581,41545539,42459286,44844803,45731150,50714388,52119702,52320471,52441025,52811289,54493573,57914481,59610309,61412083,63290799,63837536,6,70787156,71050672,71438632,72751162,75413653,76959642,82299776,82804142,8,86330733,87472850,87991393,88215621,89708798,90135922,91671884,92931994,94645214,95794629,95928781,97589886,98132477,101659947,110436602,110587467,110750000,113418969,113724414,114532808,118193661,119827719,120817779,121675348,125550255,129267451,130358669,132974642,134852330,134898638,134942076,135198169,135293746,135392405,135442967,135496247,135720295,136861914,140181623,140285895,140455782,140483073,140569927,140663378,143958421,145194841,145296337,145419284,145738761,145796012,146594475,146812752,146868942,147220034,147743475,148155154,152862259,155087696,156112073,158410955,158568795,158691051,159122335,160907269,162764958,165106451,167055484,167583135,168943170,169256352,169863642,170096772,173663306,174803358,175132920,175231120,175327531,176039673,176102849,176169366,180329893,180529468,368271,832604,1015537,1256067,1323443,6003353,6491368,17210287,19529171,24234287,28158076,28914669,29629088,29703276,29725743,29799084,29963273,30002118,30018180,30032179,30052380,30082198,30136914,30146859,30176470,30203151,30247696,30526822,30539540,31031186,31218678,31346830,31384219,31431924,31475454,31569664,31600949,32597719,32740135,32837319,32914261,32955475,33043873,33156393,33192717,33237268,33267968,33275862,33283538,33288987,33324080,34050767,39124284,40662572,41082794,41447213,41482504,46659370,49500958,50790292,50895244,51989252,54641636,55455701,63053813,66067973,66554395,69401936,70432850,71048760,71179612,72948757,73387662,73924304,75851599,75968298,7,78228949,79614333,83631972,84474315,84619538,85529420,86092457,87703971,87824125,88932115,89373710,94182992,95737264,96131961,96570589,97076231,97391550,99379114,99833192,100548126,100690600,101001727,101948145,103184180,104992219,116354921,117691276,118885677,123358862,124165875,124848782,126110414,127006295,127480957,127629525,131076467,133603778,134026813,134200818,134252331,134288828,146392214,146797167,146906173,147130357,152170514,152999502,153177757,154402278,156759762,156953203,159573222,159656679,159963628,160599366,160688403,161020081,161986906,162727667,164936154,166273926,166491834,166586826,166952319,167047645,167553920,167684511,167754854,169580182,170863413,126074,177406,212643,244576,926249,972101,1281575,1295050,1375041,1410012,1608229,1653241,1719971,1750552,1810233,2524896,3176282,4840629,4867861,4889232,4964930,5127298,9732185,10785829,16076258,19112396,19714921,20127170,25758209,27145586,27185833,27190791,27198329,27205214,27211054,27219205,32076587,32304507,32433958,36860816,38636963,38724171,45579910,46805988,48858054,49783553,50103117,50313440,52181841,54577227,54699272,63405184,68399700,70234163,71438692,75702236,81078316,81909956,82629609,82911859,84652774,86253555,86526767,86619559,88226423,89585827,90460857,92703484,93042020,93313195,95546427,98083740,98277882,100500000,101939313,102099182,102116472,103416760,104301561,109842375,113512159,117611219,119700075,120183063,122959967,126678535,126937223,131280713,131911013,132093937,136204393,140418826,140892235,141727252,142204690,143027414,143164773,143622086,145443962,147552557,149375334,149526539,150048346,150127741,150185351,150275387,153076096,153632896,154172991,154589411,155367477,155548790,155590900,155608270,155640093,155786777,156074373,157079813,157097769,157216917,157243307,157324079,157340560,157352193,157383319,157553098,157673546,157730468,157796096,157891305,157940648,157970363,158821424,822243,883458,1035965,1082156,1991215,2025112,2117269,2787051,4836375,4967724,6106395,9793586,9949038,10468642,10910733,10954071,11095093,11179303,11361188,11458744,11574175,11643003,12621393,15138607,16270766,16903414,19841122,20098825,21652357,23595751,24826812,25953116,26777558,27194965,31616066,32524948,33307911,35212220,37341975,40874320,47937993,49335553,49471615,50984822,52483714,52865800,53640096,54299722,56176209,56524409,56781719,57394923,57513500,58068253,58217498,59069554,62575831,63219500,63932587,65444456,65651089,65873543,66652437,69027137,69405473,70401120,70906651,72622117,72916427,73079858,73611212,74043287,75250000,76059082,76482004,77752721,78074686,79510101,80686529,80750000,85259313,86069241,87150768,88902009,89409080,90773919,93025319,97225938,97287918,98705869,105304559,105411388,105547847,106399693,107199114,109163660,109266967,110282497,111055288,114514001,116729605,117521395,118602035,118869513,119703106,120033126,120289568,123639212,125337303,127482727,127906835,132121384,132978958,133210746,133561579,133756754,133856622,135526578,136538120,137613038,139577976,139959062,140700277,140783766,142597366,142699308,142982329,143052690,143200560,143330953,143553715,143733653,144006276,144131451,144338703,144383082,7773697,10602634,12765485,12888008,17896418,18428555,18592436,23810690,23840909,25334004,26735223,29202170,31243815,32374419,32772935,37058454,71848655,7,73251332,76301597,77650412,78710622,80750490,82433196,84487083,87201662,95628418,95747936,99650516,102830432,103200673,105801450,108296633,115678414,117956182,118489293,119016617,119215073,119547047,121170906,122155437,124015573,124586354,128024683,132112881,134026533,134351812,134388115,135463990,135556393,136318676,136532868,136623655,136800333,136832764,137106930,137275870,137333304,137444288,137731103,139935581,140037104,140061128,140134463,140161970,140189215,140231189,140273252,1242140,1274150,1358902,1429620,1496343,1573536,1731423,2347247,2533462,2805002,3079854,5402846,6661912,6985026,10348163,11099447,11247183,11355137,14258714,17283309,18469130,18607093,19803915,21106711,23682112,24022387,24516738,26262926,26540612,31647232,43820986,44010310,44680179,45746970,46060676,46182455,48447597,49533625,49993189,50204279,50557570,50639663,50895005,51040816,51582618,52420961,57058073,57791064,59536387,61007220,63803794,64234440,64419127,65403272,65600622,66483691,69193848,71482323,77542088,78497256,79732785,80272086,81114168,81231464,83623799,85314678,85889037,85944178,85987088,86059376,87155144,88012619,88082554,90332507,90629766,90740272,96916102,98905247,99779603,99884254,100982145,101078999,101124449,102404414,102420686,102576114,102611061,102895703,102963489,1,106389556,108913769,109664185,110215916,110661713,111206593,113799878,117874802,118020721,118351165,118912340,118966038,119324702,119740455,120033075,120414301,122206248,122547783,124568115,124628732,124703257,125415484,127381162,129424399,129560921,129953956,130155095,130398432,131148389,131527087,131576354,131646679,131703031,132255674,132473709,132758193,132948327,133208017,133271394,133358343,133494541,133524257,133780584,133798314,133890187,133921960,133945476,133970436,134522358,134542840,134569923,134605652,134624202,134671722,134805497,134821706,134843110,134864692,135228872,135286136,1046330,1073133,1084945,1097801,1145475,1314201,1360789,1420870,1630409,1671922,1725554,1791172,2110608,2169417,2345820,2398645,2421746,2500000,3138149,3181099,3195792,3326514,4585673,5478047,6172608,6660871,6904095,6997647,7066707,7423119,8571972,13940641,14044611,15051532,15538944,17522341,17674080,17749851,18167459,19095072,20285349,20647154,22170747,22319437,22410875,22603281,22807340,23383527,24474799,25961581,26500000,27033864,27828181,29994247,30301111,30562133,31347662,31965377,33581788,36354501,43290014,43623296,44498112,44910465,45125193,45333506,45627772,45750995,48858272,56704521,58053066,58429426,59043273,60039358,60171241,62977130,82078232,84987972,87366026,88692580,89236691,91513930,95927232,98396752,99693100,100790166,103412614,103939335,105337542,106393532,108798728,109368791,110661709,112337733,112648918,113412129,113750000,114336964,114529816,115876391,117156425,117252311,117735373,117750000,118797898,118960362,119006993,120246625,120399898,120503595,120668453,120738951,120803151,121446579,123117107,123997873,124809653,125230282,126375435,127826870,131065565,131285537,132317871,132735993,133208261,133239910,133409544,133786589,133841501,133899580,134039451,134110291,134137254,134336538,134452384,1912416,4788533,5023272,7483672,9283725,11544498,13607115,14025892,14186787,18885877,19985403,20412882,20595791,21481735,21985220,24606222,24750000,29136803,30214029,34066496,34747961,39372788,39868421,40253531,40549466,43730468,46218811,47208040,47302799,50344020,51394178,51469979,52607567,54305222,56354739,56373797,57529059,58275917,60799842,66908041,70057189,70951949,72666508,76242583,77107933,79626164,79995699,80676450,81265011,83829330,86949170,95824932,96388132,97165422,98902373,98938712,100127517,101742249,101835169,102220219,102413277,102658651,106235732,106761566,107047079,107432831,109955443,110232560,111822420,112513700,112727583,113318293,113517185,116683149,116781477,117696492,117902973,118515876,124381507,124584052,124734392,124911851,125241513,125776730,126196463,126506403,126811544,127316993,127755687,128353687,128750493,128905635,129388312,129474729,129691800,130349045,130483438,130509471,130543420,130611057,130711656,131519816,26232225,27389487,28811885,28952341,29208838,30275144,30633108,31318629,34167576,34590946,35603039,36145948,37341681,37659658,39644368,47292563,52317897,52672264,53592326,56612971,59109123,63283719,66702593,68456904,69579732,71147084,72173959,74733096,76122448,78068114,78718711,80127047,83351664,85658325,87121569,88813412,90538851,92677245,93969295,96377189,100839335,101366425,101844680,105236573,107316334,107665107,107945798,108349600,110970232,111073562,111118824,111151995,111595458,111865816,111943427,112027798,112154337,112351306,112645553,112724781,112762784,113651130,113782454,113821717,113867101,113982201,114018249,22031905,22958983,23009542,24588263,25744002,28304639,29465960,30064279,36042919,36710982,37122395,37746904,38534690,41136711,44425877,46190617,47213182,49024691,51803956,5,55654794,56064449,56588297,57401428,57652137,59112994,59263583,59406703,59934195,60021728,60173730,62200132,66035943,69724847,69770566,75888897,75912258,78804785,79994522,85066246,87528648,91858665,91950366,93323704,93474987,94216636,95411961,95575063,96569357,96754811,97164257,98329219,100122857,100190760,100223609,100360275,100418137,100601395,100781966,101062474,101164606,101241755,103621764,103871488,103966948,103998485,104112351,104174537,105306532,105387983,105741522,105814646,105925849,105987405,106292136,106360585,23105710,24425189,25169606,25460312,25688940,26013368,26373444,28270888,28441781,30109455,30247457,30391938,31390107,31781541,35570175,36091387,45152320,51838946,52080271,53250620,5,59306912,59907948,68068338,71130831,71250000,71446912,71522374,71630688,75897895,75935017,77361884,77390274,77890313,78483130,78647886,79346258,80524340,80713731,80798714,80910871,81112170,82539596,86600525,86799817,89640951,90197016,90260269,90737680,90919936,95997024,96304537,96782129,98649747,100318876,704748,772827,833701,969878,1038522,1079732,1119618,1138521,1205151,1248004,1291520,2807245,2829440,3160438,3207043,6008914,6473019,7294333,7508195,7932588,8559699,10181443,13921251,14734695,18338665,19992207,20267267,20537310,21295435,21675830,21738982,22454917,22833843,22980314,24174540,27028511,27981317,28322357,31487570,31619339,34262015,34444217,35143302,47866623,49704990,51137136,51609792,53157686,53915021,54070720,54100420,54226696,54352161,54424147,54782072,58346292,60626621,63542181,69821860,69880615,70017281,70930814,72788677,74826477,75224479,78180991,78833268,81218151,82398912,84877408,84968654,85085771,85272498,85427362,85469903,85708069,86235385,86775653,86913222,87023701,87596094,87636386,87785252,88765240,88822254,56118,99718,185074,3236111,6288188,6399852,6557145,6619928,6738152,6816985,7983853,8809193,10367386,11084211,11401795,11789450,12350907,14855147,21220098,21296546,28642457,29391012,29508119,29930077,30312286,30724599,30838346,30918978,31622177,32363500,33412810,33592248,41216446,41327203,41382386,41794306,42011805,42139650,43210550,43222101,43231083,47590173,49189648,50332864,51826224,67623573,70178789,73300920,73395397,73465484,73500000,74690422,74971839,75307036,1227142,4443968,5186242,5373906,6274373,6403958,6577203,11138306,11679060,11740920,12244146,13252328,13377384,13446351,13631583,13858531,14047325,21184015,21273517,23018948,24009411,26765639,29274507,29993044,30056355,30327441,30810651,31074782,32579621,35113580,37789958,40358900,41609315,42394359,43529402,48122374,49900777,52204578,53170704,55126821,57151662,61568754,63334216,64191682,65218488,68359953,71756904,72191357,72220224,72282226,72454824,72928131,73090543,73509219,74099661,74170209,74224599,74564350,74702741,74740247,74775196,74825846,1474704,1497218,2653531,3226607,3285243,3819585,4558078,4710817,5241402,5286043,5395085,6541339,6610885,7676814,8668812,8772314,9981749,12967816,13059698,13226921,13258417,13270560,13477751,13674260,14220634,14237687,14305342,14661834,14951015,15523112,15556297,15694732,1,17619054,17652066,17889349,24423622,32976367,33909840,34275785,34395409,34707620,34788904,35405247,35557522,35626661,35633692,35902356,36424363,37407223,38804118,38867183,38954782,39281136,39955487,40015370,40087777,40534269,40716021,40727841,40740395,40795393,41114249,45746963,46223643,46287734,46319315,47040514,47154784,47189658,47231494,47373143,47893938,48602682,49808889,51148048,51185413,51491892,54623014,55799236,55853734,55881635,55889954,55911865,56103445,56197725,56376524,56766122,1191458,2677996,2749638,4077481,4150147,4176531,4250000,5433553,6696051,7897260,9435110,9767270,9963609,10100767,10251573,11605791,12828140,17154527,17244124,19125019,20292399,21433931,21544955,22496966,22615264,24397843,24838205,30069404,30162404,31350798,34493470,35521120,35755595,36805110,36973283,38691826,40060637,41250888,41976510,42177735,42221629,42511676,42872150,42920293,43355362,43797801,44119177,46013197,46206187,46435114,46582630,47531918,51023113,52525139,52652427,54012009,55596709,56093485,56658100,57585674,57613202,57846531,58975902,59781492,59882169,60524832,60589086,60810302,60873751,60897140,61166025,61276698,61444304,61448458,61562769,62200156,62280313,17806676,21291293,25655978,28900367,30233256,31336240,33272595,34963174,35073651,35320829,36181060,36261213,38954112,44613517,44951317,45572480,45590156,45671932,45877481,45990213,46111439,46142936,46342078,46944323,19952175,20025100,20100895,20127081,20799165,20907656,21152789,21316927,21742223,23489864,23531952,23678150,24008575,24070525,24895198,25154977,26164424,26340272,30901857,31084981,32258813,33965811,35823649,35985550,36050664,37569113,37653148,37708691,38372550,38720429,42068825,42137343,42215733,42539279,42589274,42682234,43057958,43196817,43394596,43442822,43884418,46406249,46983820,47108721,47263547,47524256,47683776,47833342,48083536,48946541],[35,37,21,21,35,49,29,26,24,50,46,34,29,33,41,24,52,32,33,40,63,40,37,28,47,23,38,21,31,49,94,93,33,33,75,37,36,245,21,59,79,54,32,43,65,21,85,80,39,58,21,68,36,38,64,22,120,130,29,88,112,88,261,88,38,240,31,69,58,105,128,37,82,48,25,101,103,22,68,70,30,69,86,76,48,61,47,87,35,50,71,72,100,135,55,97,32,79,182,64,34,69,37,25,23,31,57,73,31,63,85,72,63,73,31,58,32,31,26,42,168,150,400,64,21,67,379,22,39,649,144,73,70,48,132,65,307,43,170,185,104,95,45,78,235,105,21,32,38,101,129,44,141,216,336,28,146,101,54,312,67,35,47,79,49,135,42,28,76,21,104,96,120,26,80,124,71,62,49,145,282,43,64,37,50,191,57,38,60,61,29,41,39,75,23,24,46,72,50,39,57,82,46,82,51,70,38,53,49,38,38,43,44,30,60,26,40,69,41,92,81,169,128,35,30,64,32,74,29,195,93,64,84,42,55,50,33,86,61,57,89,86,88,50,28,25,48,22,75,128,242,30,186,105,117,30,57,21,21,55,53,35,23,116,33,37,76,22,126,68,37,76,219,92,33,62,31,85,33,33,31,121,85,52,27,132,56,69,36,143,36,98,311,23,111,98,70,60,48,73,91,23,90,159,93,98,114,49,110,103,74,28,28,41,69,27,40,30,189,85,24,29,66,24,77,61,131,49,23,48,31,79,28,103,114,51,95,68,135,59,46,52,44,91,134,176,126,90,52,43,49,23,37,57,77,29,61,71,97,38,75,85,30,103,29,89,95,53,69,30,21,36,22,69,35,75,24,67,106,25,96,27,64,38,34,42,42,96,136,56,36,39,28,28,22,83,67,27,53,70,65,22,65,86,30,22,42,34,25,172,211,125,144,33,254,54,46,66,39,309,276,39,191,53,85,52,59,78,117,81,180,192,34,117,84,383,65,485,36,93,50,151,91,94,188,78,46,223,138,26,31,137,162,202,120,56,110,44,46,130,34,152,65,24,34,32,36,38,22,35,56,25,41,71,113,45,64,63,45,56,21,26,76,176,35,42,49,87,55,28,125,168,238,123,49,43,100,35,109,195,71,24,77,99,53,137,77,54,70,42,56,59,59,81,43,84,46,72,44,81,45,70,31,44,48,37,83,46,64,43,130,34,53,78,39,109,21,48,39,73,70,146,50,54,23,28,33,84,67,56,59,43,43,92,103,90,39,36,24,58,37,132,44,46,28,35,57,44,121,106,50,65,60,133,48,31,56,56,24,24,35,46,77,41,63,34,40,26,61,26,55,41,50,43,22,46,45,59,27,23,105,37,21,88,51,51,60,53,105,96,105,39,72,28,69,46,78,57,46,66,26,38,71,46,41,28,30,77,69,28,24,94,100,46,59,31,31,64,84,25,27,71,55,81,92,50,36,71,92,21,36,36,37,27,32,21,22,70,55,36,67,59,31,41,89,44,118,97,141,25,23,198,77,64,30,37,76,37,58,29,21,27,44,82,68,54,58,33,32,23,51,24,23,49,21,60,46,59,41,21,22,82,25,35,54,36,24,49,34,60,34,38,40,22,29,43,29,88,46,43,106,30,58,152,74,57,79,147,133,107,82,68,49,93,97,63,83,39,52,29,59,33,85,49,50,32,23,21,36,52,30,40,93,31,42,101,85,80,415,109,52,47,29,55,22,53,94,97,169,28,119,51,27,103,29,31,271,40,38,22,40,230,131,181,52,315,95,148,97,211,257,120,286,25,57,52,57,43,33,45,148,102,28,48,59,143,45,44,37,33,29,26,24,87,30,63,32,24,76,56,82,73,23,60,36,23,32,37,80,25,45,24,49,42,124,34,25,22,35,58,172,31,34,31,107,51,30,29,37,43,33,22,61,23,95,22,48,37,50,64,54,32,34,25,26,25,26,21,49,40,27,56,120,26,35,35,41,54,30,47,27,34,21,63,72,48,161,58,48,44,32,63,45,43,67,66,28,26,30,63,21,51,24,67,29,58,30,58,25,37,37,31,39,60,24,79,29,138,62,81,45,58,43,48,27,22,28,29,26,119,42,31,86,23,26,74,63,31,30,41,34,53,49,41,21,76,75,122,31,28,93,81,21,21,56,37,27,39,52,89,25,50,56,41,51,63,170,23,53,24,107,44,110,104,50,68,43,65,142,33,80,95,29,98,53,93,162,39,24,109,38,70,79,22,105,49,110,34,32,125,28,23,22,42,85,44,40,28,61,32,21,28,49,21,25,24,28,42,264,31,56,43,77,47,28,80,32,46,54,40,28,62,32,46,95,32,46,55,28,55,81,29,52,37,81,54,93,39,39,112,32,32,26,74,167,53,211,104,45,101,36,108,68,83,203,60,22,27,22,78,69,65,109,43,39,42,64,39,21,60,41,111,55,33,33,90,69,25,77,59,43,32,54,40,86,24,26,106,27,67,34,127,52,63,80,29,73,40,51,47,44,69,37,67,42,58,56,51,38,63,90,36,25,45,28,146,160,98,42,52,103,70,64,32,72,103,33,51,71,22,37,43,51,61,42,39,63,79,93,71,42,40,63,21,72,86,23,36,39,21,33,79,23,36,34,22,32,32,27,55,30,25,22,21,24,23,38,32,73,35,177,131,75,22,87,53,75,41,33,33,25,69,75,156,83,41,336,56,52,31,47,140,52,31,172,69,29,53,37,29,26,46,46,21,56,76,22,116,51,27,47,97,71,132,28,145,31,30,35,27,28,38,41,21,44,56,91,134,24,21,27,60,51,60,35,130,72,51,78,109,43,31,23,31,45,22,35,21,42,54,39,119,112,25,67,48,113,24,201,56,31,52,23,56,44,95,117,21,36,42,21,112,39,49,22,51,100,111,195,117,97,189,112,70,30,103,54,33,44,45,24,40,33,66,95,64,22,25,77,25,24,39,32,26,25,44,25,23,36,54,89,158,52,21,26,28,67,64,30,87,49,38,40,33,22,23,22,57,29,36,28,103,52,31,28,38,61,75,164,23,44,38,63,39,40,56,62,64,21,22,35,22,68,22,23,24,29,46,57,25,47,68,28,32,21,23,34,24,48,50,40,24,23,21,31,61,62,38,54,33,42,78,100,142,128,32,56,137,29,41,40,65,38,267,34,40,75,31,27,31,28,21,38,54,82,40,100,40,24,43,56,47,71,26,30,27,55,42,38,37,23,32,83,120,51,45,88,175,23,50,294,68,121,83,53,34,75,135,34,179,74,70,112,190,30,255,143,160,40,79,25,63,40,145,49,65,35,66,22,101,51,59,115,327,48,58,76,208,26,45,145,55,115,30,61,167,89,146,28,54,26,25,23,221,43,292,225,76,301,189,130,42,260,34,94,25,25,23,65,28,26,30,39,176,44,26,28,70,119,24,30,56,44,31,26,27,40,33,60,64,154,41,39,53,33,97,40,56,46,38,21,27,25,29,62,30,72,147,746,68,45,426,588,288,96,209,78,340,65,199,24,27,116,101,152,170,94,230,159,124,21,33,45,32,108,61,49,23,51,74,125,24,178,163,36,99,85,83,109,42,62,78,94,67,111,50,211,267,384,48,32,34,94,39,42,42,43,22,40,26,33,36,22,27,44,70,23,27,36,280,37,52,25,72,276,117,35,34,51,34,125,87,292,24,142,234,53,61,237,79,121,29,56,21,199,122,108,158,67,191,75,55,92,64,141,92,240,178,46,48,350,198,53,33,38,83,32,39,34,36,31,24,54,43,40,57,24,30,87,25,29,26,41,40,525,50,104,265,712,128,125,41,169,98,170,101,166,202,94,200,31,63,46,41,95,158,48,147,31,22,23,29,28,34,24,279,177,70,33,57,33,539,54,320,172,97,51,180,87,128,32,49,54,120,174,107,27,52,36,21,22,95,183,67,39,66,80,88,22,42,73,48,62,80,341,103,163,127,23,63,195,30,71,211,105,55,21,28,88,73,44,23,31,208,74,64,428,238,98,414,26,67,221,91,98,54,42,26,31,84,44,67,24,104,104,128,78,47,22,26,43,34,89,659,36,52,375,49,293,126,23,30,81,49,57,22,42,35,31,32,22,41,51,30,88,54,36,41,21,43,63,51,63,40,28,32,113,110,31,40,47,31,69,45,60,24,32,30,25,99,90,99,24,65,79,26,35,94,71,142,601,125,46,46,40,23,36,51,21,46,688,338,115,200,27,48,91,91,193,163,83,21,198,538,207,39,105,108,149,46,309,25,26,22,30,31,30,21,133,121,24,21,25,79,54,71,61,45,59,22,24,44,41,55,34,61,22,24,56,30,93,202,29,178,85,54,36,28,22,21,21,25,88,58,27,23,161,89,44,21,22,570,1299,228,21,346,93,186,51,52,26,27,205,75,65,63,33,33,56,25,33,56,81,21,30,29,41,22,57,21,35,230,48,45,83,287,194,129,34,23,53,21,98,391,80,69,23,30,102,246,253,35,83,41,437,584,126,45,45,314,42,169,21,21,35,41,26,43,29,34,47,33,25,22,21,30,25,48,53,74,22,35,32,41,33,24,58,52,31,23,59,27,41,76,100,29,29,30,79,38,54,28,40,61,57,32,49,95,39,24,38,43,27,97,343,139,47,61,24,61,115,61,25,194,550,393,103,78,100,120,36,66,36,21,63,25,25,65,41,25,178,85,34,197,21,35,37,37,127,394,26,28,41,216,26,56,31,40,21,30,79,158,25,96,26,89,22,26,21,23,26,32,52,98,55,52,27,23,41,61,26,55,28,58,99,37,61,64,62,42,102,26,47,191,26,53,26,22,41,346,162,115,30,107,49,21,38,29,69,151,21,83,69,90,27,75,96,66,42,138,27,215,54,70,36,88,171,44,40,32,24,26,49,43,47,26,42,38,23,26,52,94,215,91,106,26,291,24,59,38,61,342,49,111,37,145,38,25,22,28,67,70,21,81,30,90,46,32,41,78,75,84,35,150,71,33,52,31,50,201,122,208,32,74,147,230,355,28,36,55,40,49,41,78,290,80,49,25,245,21,25,122,51,78,27,72,203,314,108,93,25,49,25,21,52],["19922","20362","8728","5427","16819","13768","18140","39032","10166","54891","50422","18760","25670","31440","40230","6552","28199","18196","12629","21246","35934","18176","35460","34766","31783","10232","19271","15763","11110","197363","492315","583114","80645","19797","55756","30196","100209","600671","28235","135870","112143","51025","227158","58519","222444","24963","281261","254199","105375","147226","70084","66824","70011","42151","116926","44009","146679","799853","260301","113190","350550","565835","491201","282985","25990","393492","129844","151452","57755","258876","358413","304686","1124545","284345","81378","605933","824633","61901","750033","1069959","350546","534810","1201686","1914088","502332","356682","187633","792456","227702","229210","514763","1314834","656326","669685","903247","601950","259012","354701","2019217","2531076","975349","294658","147875","31827","50619","66056","105284","182446","113745","75989","321370","91238","189162","170818","62903","113265","86546","133080","61354","60231","298903","244686","532210","71157","16166","303605","395094","14485","27727","821813","575571","74780","59805","38981","273529","126803","790031","183246","671891","478194","586457","469534","141996","198902","900013","718417","53921","94415","68461","312230","215483","165931","392815","4053965","3179555","54577","142723","153050","75039","632218","205731","281919","110398","212046","436353","444420","227342","239113","391638","288296","651959","406620","264375","35125","68156","215373","212535","289982","308425","287801","1089832","533192","364913","45897","198784","754302","52739","143252","130832","80593","105027","110284","95196","194114","21931","27897","55973","76628","52558","45134","247638","460390","100382","289752","87833","51361","81886","216329","848996","647615","285245","335686","260280","122291","77384","24922","47077","112835","1255471","318151","279343","2289356","681702","187118","213017","2263657","199062","1096568","125158","154288","1477108","679703","1542903","793491","946912","83741","295488","801632","2109399","519102","789435","871404","2823031","25836","11162","132217","329136","113770","139941","201639","275703","55136","1642018","250661","656483","280233","136734","47341","32883","152345","61210","37239","29109","636196","411788","1934856","387002","152315","144728","66640","80540","88181","2246337","2630439","85566","164223","334208","94940","52751","30385","167668","570547","369633","192166","122800","304008","341157","1054787","121534","3227327","542062","749500","2062530","243411","872454","507318","611706","392446","824085","1302041","691496","92587","427395","496744","385734","916320","691094","209131","453810","972915","476245","198190","921451","1139801","1227877","1520548","66962","17478","1555615","197865","122454","468336","516710","641514","1090860","115479","240480","50046","38051","25758","9127","1221082","103824","561014","539250","732533","254199","828854","169446","61376","150073","84448","100246","101510","212551","114009","116408","83352","87820","45898","81834","27730","75678","36025","47167","16864","34802","40441","70913","110518","1789134","924020","314715","1545839","29573","146590","310296","175997","136848","157112","131979","77397","17713","701188","289620","1294228","487198","688358","498088","225806","740072","322983","282266","306821","235168","155363","237227","211500","495299","402797","265249","507667","124739","366928","139584","923614","884776","88017","989129","1132188","1107891","185524","2151038","2030771","137324","564952","221700","900603","1340749","1838456","942862","408947","708997","523386","3055214","130603","197910","794776","510141","2212517","890128","109565","507209","234525","201986","117077","72407","293098","507857","191673","1066308","557157","237638","249470","303304","2067684","705308","1273049","101380","528556","550827","506714","2572118","1915468","1308801","354176","269005","596719","512022","95149","519711","820683","564576","705385","273789","187197","470156","343996","395261","358510","104076","175674","193155","21621","13816","30268","21601","13921","13079","24053","23934","76582","43453","79820","74286","146828","475961","114750","97016","86141","20840","42666","90215","288072","80939","44101","40166","85380","55001","27584","766308","1896629","1234851","789595","1564686","444947","1133103","365490","1074917","3687982","1087225","976585","3394272","998859","257474","1789834","1663182","602554","329320","129890","108212","175184","267410","210889","2269150","1831597","2319607","1755751","230655","376924","1218764","491012","276429","535817","332252","74323","392878","323880","643987","159489","535082","274108","1063925","1522486","285575","2528638","282464","771856","648130","414342","1304846","684887","618610","1999496","219975","319070","633409","1143139","2570391","1169970","1456923","431149","290065","535037","885723","539888","83291","57738","76300","425838","162183","3544080","164806","578831","313445","382971","1528257","973883","617274","466709","399786","398826","2366287","1361041","224211","219381","647796","281013","163922","40803","24302","36155","57747","31382","44763","15289","12333","15116","42833","14234","22780","40507","29928","27253","10521","34202","22262","72885","24311","42944","384028","78351","18004","321092","117210","208983","54287","247399","903803","263654","490761","69905","133251","510430","690315","206656","565368","437679","518027","863124","140111","403037","427163","241746","256066","654058","1086603","2149332","916540","718227","759895","2139494","2265313","434656","338336","43763","47212","42602","417097","184820","459285","1041372","885700","412893","1388592","188688","119778","369764","243572","46240","384657","412082","782997","339673","161638","387658","47092","387160","36984","314687","135443","549775","54141","193644","2330732","668779","56900","194054","1487452","245209","871462","1220109","530319","227145","133944","839885","519590","1098498","309089","131297","161654","418968","304509","805713","1463662","845961","988685","129467","521319","149281","1089157","74943","30138","43094","34357","255288","77336","97988","50031","51872","50072","717159","12452","41556","128666","20428","8117","55820","393163","1236163","101004","122736","32715","56936","265183","217998","54810","350000","198445","140880","1577241","696782","989052","1047544","102216","67931","430917","1082742","1856995","1725934","716988","364689","282316","259071","176661","232722","314113","457316","114557","73834","94858","32443","62721","65963","79892","54459","29480","61280","182628","179619","60480","49054","487704","132641","889423","1555602","188501","24387","387135","73463","17102","73130","23207","37573","15083","12832","19867","29049","53565","9684","24990","23743","44106","26821","11930","20760","186913","72927","36119","47418","42263","89767","26781","325357","79825","96313","76405","25126","73109","72488","35924","44294","30271","6546","2736","4621","34845","50766","100160","162571","407235","179739","31982","92153","799063","1040291","94909","459870","557800","308349","817264","1168893","486192","2540383","1030645","298759","129978","295631","438063","534473","1282157","116489","26892","1469053","1197269","499933","196255","143440","534972","551019","1258093","119698","117516","439541","2805051","1550865","194696","437825","505373","314728","1541193","330615","374127","141726","217907","424492","1226484","1515889","1329181","191275","537414","108861","806422","681773","383849","303214","474275","146071","347037","288741","421214","173664","51288","29854","64819","404693","108578","168418","46601","191790","177591","877984","1082187","192963","61914","81291","306664","99142","88711","315955","153954","392895","186153","85474","149602","83412","163098","85975","46370","130380","69739","80181","99712","65098","50891","34967","28734","10554","11000","27390","12576","37651","19432","31665","25933","42927","30319","34346","24895","119749","16936","25105","20952","74827","37607","1281948","1053088","1923899","526675","563362","347639","508208","8422","4911","4515","2760","5591","4217","7382","733519","227071","128707","360815","169006","86637","415577","441634","392754","848658","316818","143907","681840","899572","116980","114786","1899699","234162","1202373","32757","328315","831311","718152","281510","245270","1599376","272718","92384","340334","1358607","396556","191208","211556","270877","306058","145514","192186","103486","16732","54472","16903","542683","544365","1766165","1144705","239601","2058117","480341","613316","2321099","116902","280712","629976","180869","768100","57413","471205","491394","234761","27413","136646","98043","1185553","2107539","125333","149136","271850","78789","56211","89712","821857","251518","171423","95232","50821","181087","41789","17122","31618","145361","116590","17657","17685","20011","23804","12015","16183","11424","20462","157051","23511","40679","41588","87662","27177","27104","71423","41444","43052","56695","40944","31401","33645","51434","318723","1809599","127680","1138040","43585","120079","91757","175688","42486","128471","82650","118088","96636","65516","38381","53795","1852261","595652","215486","48766","256538","1356743","53517","1186105","363434","199545","414481","454825","907059","781467","1512077","2128321","352961","289453","85552","15514","931394","1498300","381025","390095","282742","479622","346093","256996","63269","118150","544883","148491","714115","360845","268241","607968","944455","197588","210843","778160","431552","377972","276491","183636","191773","284311","160636","284466","431007","54401","633482","422332","1269969","317270","1434488","806570","61679","2342281","809438","240802","1152008","409079","1364321","650066","297378","44861","346288","721469","106232","98684","728812","797817","176343","101941","413473","282214","3456752","2215343","770597","446626","267020","509492","329631","255749","347416","292984","701797","90247","345426","855026","77607","350580","193985","99441","599362","288119","840756","1332052","380105","740647","83251","75673","101096","142264","70129","147604","130177","10534","40484","92868","25174","25214","43882","507539","1753925","1156896","122194","326613","531066","163633","1372994","28241","1492481","985222","1243610","1427078","1129655","209707","30112","206449","150902","1026063","1317762","431327","210621","491466","266963","504029","440703","128417","84633","44515","330431","368678","2260552","546632","160980","1063313","532296","256198","198241","329729","1623442","983401","53581","336353","206883","44496","81376","76093","35812","22892","91823","66446","34244","90528","126114","32200","232633","67754","36499","110661","158556","38352","100494","23254","73115","27278","24280","41176","41840","20925","31644","84479","23370","66328","66222","145630","577228","185990","249779","79853","289272","218146","322377","857947","750851","146652","107150","284523","283308","404448","136852","447750","764557","153562","215140","491753","110659","276583","147231","70985","189033","430178","105373","163705","121218","169681","50266","458581","194094","64897","81794","37588","145592","84324","250933","2558072","729842","1745105","214677","53793","430013","170342","153271","197039","882794","2709780","249440","544972","954944","221797","538491","114167","100714","1123798","1688968","574054","54557","41804","71968","655143","182923","69294","332506","296564","109668","473321","48876","179496","102943","692888","94820","43803","91664","11027","76477","30668","5449","66805","13363","301009","2520966","749389","550961","445257","544276","970584","1030717","145701","326944","22020","40867","21148","255473","120653","68489","389974","340415","175160","60298","73473","415483","381161","191101","133564","54345","200826","168664","749740","44480","20048","59276","41344","255673","217737","284038","77640","175184","63079","49857","135959","29452","51338","10166","20420","30185","23292","24625","48092","9932","20964","14628","17359","14286","10714","15915","12060","11718","36320","45167","19502","26391","11586","11281","47433","26126","43806","32459","77689","41121","52975","41171","96715","24595","78737","5417","22767","76605","39695","37754","14138","115461","335672","891940","256905","27198","242575","92858","68181","192307","205017","172046","101685","99047","445822","67138","151451","34988","83171","324416","143275","67178","1522325","147540","90921","192186","162146","575583","765861","1486266","189611","356271","127040","1489362","304240","260610","782957","160024","90546","92280","6801522","63504","199778","204127","125192","68954","293952","106923","701581","1305824","455721","179425","235099","616924","131619","429208","3018893","1766540","616025","692579","75933","1566838","143710","2453957","1296068","284518","445200","399066","350909","516069","275615","568597","548095","587732","309217","141697","174672","426129","192551","244631","261743","78959","205968","14421","40123","161208","45680","222654","37513","103302","134315","69788","63945","618066","86211","880117","268016","165689","461159","1447850","627940","218713","1031128","277475","471405","30991","76876","23747","53371","57867","39450","22499","26180","198340","115010","24407","38532","131839","213954","33724","294497","290359","417590","160277","258996","235402","427256","181401","884463","232176","1675109","142489","504494","385328","582178","353079","467607","495051","384288","295723","173365","119818","165786","94402","72468","112598","75240","175829","1036363","41415","18755","779058","485216","2523067","408040","766864","325196","1712689","242582","864361","320235","140053","680199","588205","1702575","2182126","324931","562713","620174","1088867","35998","491763","224465","92350","336459","192728","244677","176725","67865","225205","384848","40684","272934","324507","111209","213447","429532","146222","399317","97806","151664","206196","258898","143441","202281","150073","177239","329028","534581","417980","174934","304683","505203","289668","357890","396379","154829","56811","84692","172017","133483","133285","25432","33705","67260","100318","18755","201547","119718","518986","140039","140591","128525","254002","645594","427722","207587","353038","241494","177512","317429","394367","943129","105850","348267","917979","2055891","2002191","3227269","3385912","1753443","1122426","1566863","835898","1125493","762126","524855","587292","210216","3224217","2306426","1413243","1685525","1725074","1827527","1290416","835025","254376","116424","477183","2736572","816333","346042","276468","402579","108208","36667","45028","32864","43462","56701","46204","84088","108061","42379","44778","5588","9004","35015","15600","17077","7418","37045","34306","875356","29522","50293","606303","1154549","2168220","890735","597107","525520","504882","385150","584976","739401","2136710","3279797","1331487","1022272","1809351","197483","148280","222104","408461","233985","401427","249262","43390","96276","142647","69574","87297","121340","546099","540198","142099","44735","138896","23125","1287307","494521","3996058","2458416","108664","89900","356484","149500","505234","161359","161739","69356","184881","408557","700629","17874","67582","32490","88758","53950","181781","179828","66725","40689","70405","58639","111460","94925","29830","100839","32028","46626","81144","336556","72900","110966","58572","230758","68241","353562","425188","370003","189705","228284","90606","51338","87017","136043","359454","137104","143981","42980","390245","318874","520973","1368151","953177","240974","1012219","55636","160063","598939","389250","254081","117577","48329","73305","108025","183477","35108","111883","26948","377614","151412","163437","132500","89107","28109","55994","89151","105451","424820","2350524","197795","140950","555671","61582","476997","180410","129849","307141","414456","149746","196203","8050","22236","18110","8573","9740","9336","24297","18577","10068","39397","40612","18098","21436","11023","23458","778526","462616","820348","213585","288496","453129","934262","1016059","29538","590308","114198","274460","269310","74600","78076","62855","98387","99882","112936","418936","1416378","281206","34881","30884","130814","108404","144341","475926","463857","272166","1389374","470680","274481","266572","147196","29392","125920","60201","71706","357572","926304","2280047","2914045","429060","58454","67030","180813","577800","551675","397784","325269","164044","456418","1179597","312219","73258","116812","101175","154535","42256","237897","29195","25652","43167","36460","34299","37980","11790","64813","56529","24249","14303","13733","336987","288187","100172","55617","61732","117717","78536","13037","59192","243769","431129","315653","346396","350906","197740","59092","35911","237515","746433","116599","421683","322739","129958","26387","79917","87082","123403","44243","36514","92011","108511","51549","44305","205449","115800","44162","11079","8339","795550","1598182","1142964","548085","1342298","216576","201467","94248","69560","32452","44651","265974","82142","327567","573880","740594","86806","388123","129342","172329","388305","539123","60811","160303","43283","124112","68664","92773","42007","188298","810756","88614","321491","555139","2753465","668206","718198","62834","268567","482360","199094","447135","1711953","2114130","1302357","609314","118680","484725","872373","1778019","388254","230511","34300","1271776","1739986","1764595","856687","633719","2069423","459995","433710","28388","38446","117197","28552","116128","416408","111871","70097","50052","64983","44443","12489","34463","30662","9669","19536","51619","46614","40694","30361","21905","35835","39733","40831","103023","55259","68812","23099","70613","102519","43304","38663","62528","18574","31258","11914","206619","195992","38989","16714","67380","63667","201014","42498","32938","138196","116064","29406","31655","43215","276011","552744","817166","298934","119304","311448","75696","155246","145711","68330","5225","268334","521638","867173","99341","60989","87173","191216","38330","59254","65755","41397","180993","11197","11947","53879","30819","19682","145213","63842","31334","319972","17688","34122","29559","43618","143937","708528","112202","11999","37069","217607","5206","37823","18935","17799","8032","20705","72122","89663","72535","148896","36338","55297","20140","77480","48748","25765","71503","187754","554175","1198302","437091","322376","195747","136852","103188","692719","739268","495411","86501","411318","286609","107685","93023","746965","99383","480409","438167","20964","57956","350797","112970","52502","95339","11842","75021","1579310","305263","1189448","164942","199183","43191","122533","59297","47332","105361","264850","24821","164282","192612","171628","135252","162826","484405","184020","126026","1270054","63052","335824","134456","81557","26959","232348","602396","228384","75234","39539","26737","24359","62085","23101","19762","17583","19450","3716","6380","10695","71797","868489","2751547","2922040","1019064","287412","1102262","81370","53505","108556","135271","859930","79740","318466","29563","149359","25244","17411","21269","41889","101476","120084","23291","59910","31788","83325","31302","30713","24638","131831","104911","244876","85228","425008","77751","41739","145596","103032","61251","395197","258990","681273","175641","104688","182455","474307","1317173","28371","51850","64878","49047","54041","55216","87935","308092","115213","67113","54626","311545","49528","30387","173119","88858","121882","47669","99381","418262","373506","124549","154490","48781","150561","70677","59642","74615"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>chr<\/th>\n      <th>start<\/th>\n      <th>end<\/th>\n      <th>n<\/th>\n      <th>pmd_width<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

```r
# define plotting function
plot_region <- function(c, s, e, span = 0.1, cpg_size = 500,
                        filtered = TRUE, array = 'EPIC') {
  
  if (filtered) {
    region1 <- anno %>%
      filter(chr == c, between(start, s, e),
             cpg %in% rownames(betas_filt),
             !cpg_id %in% c('island', 'shore'),
             !grepl('promoter', genes_id)) %>%
      pull(cpg)
  } else {
    region1 <- anno %>%
      filter(chr == c, between(start, s, e),
             cpg %in% rownames(betas_filt)) %>%
      pull(cpg)
  }
  
  if (array == '450K') {
    region1 <- intersect(anno_450k$Name, region1)
  } 
  
  # subset to cpgs in region
  region1_data <- betas_filt[region1, ] %>%
    
    t() %>%
    as.data.frame() %>%
    bind_cols(Sample_Name = rownames(.), .) %>%
    
    # melt
    gather(cpg, beta, -Sample_Name) %>%
    
    # add trimester and sex
    inner_join(pDat_filt %>% select(Sample_Name, Trimester, Tissue)) %>%
    
    # add cpg coordinates
    inner_join(anno %>% select(cpg, chr, pos = start)) %>%
    
    # set a grouping variable to control alpha for placenta/troph
    mutate(alpha_group = if_else(Tissue %in% c('Trophoblasts cs', 'Villi'),
                                 1, 0.6)) %>%
    
    # remove second trimester data
    filter(Trimester != 'Second')
  
  # get pmd data
  pmd_data <- anno %>%
    filter(chr == c, between(start, s, e), 
           !is.na(pmd_id)) %>%
    select(pmd_id) %>%
    separate(pmd_id, into = c('chr', 'start', 'end'), remove = FALSE) %>%
    mutate(chr = factor(chr, levels = paste0('chr', c(1:22, 'X'))),
           start = as.numeric(start),
           end = as.numeric(end)) %>%
    arrange(chr, start) %>%
    distinct()
  
  ggplot(region1_data) +
    geom_line(stat = 'smooth', span = span, method = 'loess', se = FALSE,
              aes(x = pos, y = beta, color = Tissue, alpha = alpha_group), size = 1.25) +
    geom_rect(data = pmd_data, 
              aes(xmin = start, xmax = end), 
              ymin = 0.025, ymax = 0.05, fill = 'black') +
    geom_rect(aes(xmin = pos - cpg_size/2, xmax = pos + cpg_size/2), 
              ymin = 0, ymax = 0.025, fill = 'red') +
    facet_wrap(~Trimester, ncol = 1) +
    scale_color_manual(values= color_code_tissue[unique(region1_data$Tissue)]) +
    scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(limits = c(s, e)) +
    scale_alpha_identity() +
    labs(x = paste0(c, ':', s, '-', e))
}
```

plot regions


```r
plot_region(c = 'chr1', s = 1000000, e = 1130000, span = 0.25,filtered = T) 
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
plot_region(c = 'chr1', s = 1000000, e = 1130000, span = 0.25, filtered = F) 
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
plot_region(c = 'chr4', s = 6400000, e = 6600000, span = 0.15, filtered = T)
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```r
plot_region(c = 'chr4', s = 6400000, e = 6600000, span = 0.15, filtered = F) 
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-4.png)<!-- -->

```r
plot_region(c = 'chr7', s = 850000, e = 1025000, span = 0.15, filtered = T)
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-5.png)<!-- -->

```r
plot_region(c = 'chr7', s = 850000, e = 1025000, span = 0.15, filtered = F) 
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-6.png)<!-- -->

```r
# Figure 1d
plot_region(c = 'chr21', s = 28000000, e = 45000000, span = 0.05, filtered = T, cpg_size = 5000,
            array = '450k')
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-7.png)<!-- -->

```r
plot_region(c = 'chr21', s = 28000000, e = 45000000, span = 0.05, filtered = T, cpg_size = 5000) 
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-8.png)<!-- -->

```r
plot_region(c = 'chr21', s = 28000000, e = 45000000, span = 0.15, filtered = F, cpg_size = 5000) 
```

![](2_7_PMDs_files/figure-html/unnamed-chunk-8-9.png)<!-- -->
