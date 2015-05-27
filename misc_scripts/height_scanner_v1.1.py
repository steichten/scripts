#plant scanner 1.1
#SRE
#24-1-2014

# this will record the barcode height measurements of plant samples 
# by first getting the pot barcode (anything not starting with 'h.'), 
# then measures two height codes (starting with 'h.') for cases in 
# which the height is in between two 1cm intervals. The average
# height and time of scan are recorded in an output.

# will also accept scans for flowering as well as senescence info
# all data written to tab-delim text file for downstream analysis


from datetime import datetime
number_of_plants = 0
step=1

#startup
print "--------------------------------"
print "- PLANT BARCODE SCANNER v1.1   -"
print "--------------------------------"
print "- borevitz lab                 -"
print "- 2014                         -"
print "--------------------------------"
print "- type 'end' to exit anytime   -"
print "--------------------------------"

print '\n'

#define file UID
name=0
while name !=1:
	filename = raw_input('Input the filename record:\n')
	if filename=='':
		print 'No input provided! Try again\n'
		continue
	name=1
	
while step != 0:
#step 1 get plant ID
	if step == 1:
		plant = raw_input('Scan plant %i:\n' % (number_of_plants + 1))
		if plant=='':
			print 'No input provided! Try again\n'
			continue
		print plant
		if plant[0]=='h' and plant[1]=='.':
			print 'You scanned a height! Try again\n'
			continue
		if plant=='end' or plant=='END' or plant=='End':
			print 'Recorded height, and traits of %i plants' % (number_of_plants)
			step = 0
			continue
		step = 2
#step 2 get height 1
	if step == 2:
		h1 = raw_input('Scan first height measurement:\n')
		if h1=='':
			print 'No input provided! Try again\n'
			continue
		print h1
		if h1=='end' or plant=='END' or plant=='End':
			print 'Recorded height, and traits of %i plants' % (number_of_plants)
			step = 0
			continue
		if h1[0]=='h' and h1[1]=='.':
			h1=h1[2:]
			step = 3
		else:
			print 'Height not scanned! Please try again'
			continue	
#step3 get height 2	
	if step == 3:
		h2 = raw_input('Scan second height measurement:\n')
		if h2=='':
			print 'No input provided! Try again\n'
			continue
		print h2
		if h2=='end' or plant=='END' or plant=='End':
			print 'Recorded height, and traits of %i plants' % (number_of_plants)
			step = 0
			continue
		if h2[0]=='h':
			h2=h2[2:]
			step = 4
		else:
			print 'Height not scanned! Please try again'
			continue
#step 4 get flowering state (f.1 ==yes or f.0==no)	
	if step == 4:
		flower = raw_input('Scan flowering state:\n')
		if flower=='':
			print 'No input provided! Try again\n'
			continue
		print flower
		if flower =='end' or plant=='END' or plant=='End':
			print 'Recorded height, and traits of %i plants' % (number_of_plants)
			step = 0
			continue
		if flower[0]=='f' and flower[1]=='.':
			flower=flower[2:]
			step = 5
		else:
			print 'Flowering state not scanned! Please try again'
			continue	
	if step == 5: #step 5 get senescence state (s.1 ==yes or s.0==no)
		death = raw_input('Scan senescence state:\n')
		if death=='':
			print 'No input provided! Try again\n'
			continue
		print death
		if death =='end' or plant=='END' or plant=='End':
			print 'Recorded height, and traits of %i plants' % (number_of_plants)
			step = 0
			continue
		if death[0]=='s' and death[1]=='.':
			death=death[2:]
			step = 6
		else:
			print 'Senescence state not scanned! Please try again'
			continue		
	if step == 6: #traits collected, write them out
		h_avg = ((int(h1) + int(h2)) / float(2)) #make height avg
		time=datetime.now() #format date as year:month:day:hour:min:sec
		f_time = str(time.year) + ':' + str(time.month) + ":" + str(time.day) + ":" + str(time.hour) + ":" + str(time.minute) + ":" + str(time.second)
		print 'Plant barcodes scanned:\n'
		print '(time--plant_ID--height1--height2--height_avg--flowering_state--senescence_state)\n'
		print '%s	%s	%s	%s	%s	%s	%s\n' % (f_time,plant,h1,h2,h_avg,flower,death)
		number_of_plants = number_of_plants + 1 #add to plant count number
		#write out results, tab-delim to file YEAR_MONTH_DAY_fileUID_plant_heights.txt
		file = open('%i_%i_%i_%s_plant_heights.txt' % (time.year,time.month,time.day,filename),'a')
		file.write(str(f_time) + "\t" + str(plant) + "\t" + str(h1) + "\t" + str(h2) + "\t" + str(h_avg) + "\t" + str(flower) + "\t" + str(death) + "\n")
		file.close()
		#start again
		step = 1
			
		