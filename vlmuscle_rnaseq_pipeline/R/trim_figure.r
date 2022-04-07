

  print()
  


b <- data.frame(a)
b2 <- stack(b,keep.all = TRUE)
b2$subject <- rownames(b)
b2 <- b2 %>% 
  separate(ind, into = c("type","base"), sep = ".X")
class(b2$base) <- "numeric"
b3 <- arrange(b2, base)

gd <- b3 %>%
  group_by(base,type) %>%
  summarise(min = min(values),max = max(values),mean = mean(values) ) 
gd

fig1 <- ggplot(b3, aes(base, values, color = type))+
  geom_line(aes(group = subject)) + 
  geom_line(data = b3)

fig2 <- ggplot(b3, aes(x = base, y = values, color=type, fill = type)) +
 geom_point()+
   geom_smooth(method="loess") +
  labs(
    title = "Quality of reads",
    x = NULL,
    y = "value",
    color = NULL
  )

fig3 <- ggplot(gd,aes(base,mean, color = type, fill = type)) + 
  geom_point()+ geom_polygon(alpha = 0.5)

fig4 <- ggplot(gd, aes(x=base, y=mean, color = type)) +
  geom_line(aes(x=base, y=mean, color=type)) +
  geom_ribbon(aes(ymin=min,ymax=max,fill=type),color="grey70",alpha=0.4)




  