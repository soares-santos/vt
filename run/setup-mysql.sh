
## Usage: source setup-mysql.sh
## This script will create and initialize the mysql database for your jobs ##

## define some stuff

dir=`pwd`
DBDIR=$dir/jobsdb
TMPDIR=$dir/tmp
DBNAME=production_run

run_mysql_command() {
    echo "$COMMAND" > tmp.txt
    mysql < tmp.txt 2>err.txt
    n=`cat err.txt | wc -l`
    rm err.txt tmp.txt
}

write_cnf_file() {	
    f="$HOME/.my.cnf"
    if [ -f "$f" ] ; then mv $f $f.bck ; touch $f ; fi
    echo "
[mysqld] 
datadir=$DBDIR
socket=/tmp/mysql.sock
symbolic-links=0

[mysqld_safe]
log-error=$TMPDIR/mysqld.log
pid-file=$TMPDIR/mysqld.pid
" > $f 
}

## setup ups product

setup mysql

## get to work

COMMAND="use $DBNAME;"
run_mysql_command
if [ "$n" != "0" ] ; 
then 
    echo $n
    echo "db not found. try to create db."
    COMMAND="create database $DBNAME;"
    run_mysql_command
    echo $n
    if [ "$n" != "0" ] ; 
    then 
	echo "create db failed. restart server and try again."
	mysqld_safe & 
	sleep 5
	run_mysql_command
	if [ "$n" != "0" ];
	then
	    echo "failed again. try to config mysql from scratch." 
	    mkdir $DBDIR 	
	    mysql_install_db --datadir=$DBDIR 
	    write_cnf_file
	    mysqld_safe & 
	    echo "create user $USER@localhost;" > tmp.txt
	    echo "create database $DBNAME;" >> tmp.txt
	    echo "grant all on *.* to $USER@localhost;" >> tmp.txt
	    mysql -u root < tmp.txt
	    rm tmp.txt
	    COMMAND="use $DBNAME;"
	    if [ "$n" != 0 ] ;
	    then 
	        echo "Error: setup-mysql.sh failed"
	    fi
	else
	    echo "Success: mysql configured and db created"
	fi
    else
	echo "Success: mysql db created"
    fi
else
    echo "Success: mysql db is set"
fi



	
       