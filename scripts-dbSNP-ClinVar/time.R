out.files = grep('.out', list.files('.', 'mhcut-dbsnp-clinvar-deletion'), value=TRUE)

times = sapply(out.files, function(out.file){
    log = scan(out.file, '', sep='\n')
    start = grep("Start", log, value=TRUE)
    start = gsub('.* (.*:.*:[0-9]*) .*','\\1',start)
    start <- strptime(start, "%H:%M:%S")
    finish = grep("Finish", log, value=TRUE)
    finish = gsub('.* (.*:.*:[0-9]*) .*','\\1',finish)
    finish <- strptime(finish, "%H:%M:%S")
    as.numeric(finish-start)
})

summary(times)
sum(times)/60
