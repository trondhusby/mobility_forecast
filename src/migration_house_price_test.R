## title: check correlation between house price and migration
## date: 21.03.2020
## author: trond

## house keeping
library(ggplot2)
library(cbsodataR)
library(data.table)


## read data
house_price <- data.table(cbs_get_data('83906NED')
                          )[grepl('MM', Perioden),
                            .(Perioden,
                              prices = as.numeric(GemiddeldeVerkoopprijs_7),
                              sales = as.numeric(AantalVerkochteWoningen_4)
                              )]

##cbs_get_meta('84547NED')
migration1 <- data.table(cbs_get_data('37943ned')
                         )[grepl('MM', Perioden),
                           .(Perioden,
                             moves = as.numeric(TussenGemeentenVerhuisdePersonen_9) + as.numeric(BinnenGemeentenVerhuisdePersonen_10)
                             )]

migration2 <- data.table(cbs_get_data('84547NED',
                                     RegioS="NL01  ", ## country-level
                                     Geslacht="T001038", ## gender-total
                                     Leeftijd31December="10000") ## age-total
                         )[!Perioden %in% migration1$Perioden & !grepl('JJ', Perioden),
                           .(Perioden,
                             moves = as.numeric(BinnenGemeentenVerhuisdePersonen_1) +
                                 as.numeric(GevestigdInDeGemeente_2)
                             )
                           ]
migration <- rbindlist(list(migration1, migration2))


plot_dt <- melt(merge(house_price, migration, by = 'Perioden'),
      id.vars = 'Perioden'
      )[,
        ':='(
            idx = value/value[Perioden == '2010MM01'],
            idx_rollmean = frollmean(value/value[Perioden == '2010MM01'], 12)
        ),
        by = variable
        ]


ggplot(plot_dt,
      aes(Perioden, idx_rollmean, group = variable)) +
    geom_line(aes(col = variable)) +
    #facet_wrap(~variable, scales = 'free', ncol = 1) + 
    theme_bw()
