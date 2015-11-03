	parameter (kwbank=69000,kwwork=5200)
        common /gcbank/ nzebra,gversn,zversn,ixstor,ixdiv,ixcons,
     &                  fendq(16),lmain,lr1,ws(kwbank)
        dimension iq(2),q(2),lq(8000),iws(2)
        equivalence (q(1),iq(1),lq(9)),(lq(1),lmain),(iws(1),ws(1))
        equivalence (jcg,jgstat)
        common /gclink/jdigi ,jdraw ,jhead ,jhits ,jkine ,jmate ,jpart ,
     &                 jrotm ,jrung ,jset  ,jstak ,jgstat,jtmed ,jtrack,
     &                 jvertx,jvolum,jxyz  ,jgpar ,jgpar2,jsklt

