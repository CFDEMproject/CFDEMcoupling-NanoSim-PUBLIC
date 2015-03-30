
/////////////////////////////////////////////////////////////////////////////
                          // MEMBER functions
/////////////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   Read all necessary global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::read()
{
    //Allocate or Set Data arrays
    allocateMe();

    //TODO
}

/* ----------------------------------------------------------------------
   Write all necessary global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::write()
{

    //TODO
    //Write particle information to file
//    OperationProperties op(OPERATION_OUTPUT,false,false,false);

    //output().write_screen_one("DataStorage is now writing some output");
//    data().write(op);
}


/* ----------------------------------------------------------------------
   Scatter global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::scatter()
{
    if(comm().nprocs() > 1)
        error().throw_error_one(FLERR,"TODO: implement daat scattering if required");

   //TODO: in CustomValueTracker
   //make send buffer
   //pack buffer
   //unpack buffer
}

/* ----------------------------------------------------------------------
   Scatter global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::parallelize()
{
    if(comm().nprocs() > 1)
        error().throw_error_one(FLERR,"TODO: delete per-particle data if not in my subdomain");

}

/* ----------------------------------------------------------------------
   Intitialization phase - done before tun
------------------------------------------------------------------------- */

void DataStorage::init()
{

    return;
}

void DataStorage::allocateMe() const
{

    printf("DataStorage: will allocate memory for %d operations in operationContainer(). \n", operationContainer().operationCount());

    output().write_screen_one("\nDataStorage allocation completed.");

}

/* * * * * * * * * * Use these functions only for getting data from the RMA...when possible * * * * * */

void DataStorage::getCellCoord(int* buf, int cell)
{
   int ncells=selectorContainer().ProcCell(comm().me());
   int buf1;
   
     for (int n=0; n<ncells;n++)
      {
      MPI_Win_lock(MPI_LOCK_SHARED,comm().me(),0,*cellwin_);
        MPI_Get(&buf1,1,MPI_INT,comm().me(),n,1,MPI_INT,*cellwin_);
      MPI_Win_unlock(comm().me(),*cellwin_);
        if (buf1==cell)
         {
           buf[0]=n;
           buf[1]=comm().me();
           return;
         }
      } 
      
      for (int i=0;i<6;i++)
      { 
        
        if (*input().getNp(i)!=-1)
        {
           ncells=selectorContainer().ProcCell(*input().getNp(i));
          
         for (int n=0; n<ncells;n++)
         {
           MPI_Win_lock(MPI_LOCK_SHARED,*input().getNp(i),0,*cellwin_);
           MPI_Get(&buf1,1,MPI_INT,*input().getNp(i),n,1,MPI_INT,*cellwin_);
           MPI_Win_unlock(*input().getNp(i),*cellwin_);
         if (buf1==cell)
         {
           buf[0]=n;
           buf[1]=*input().getNp(i);
           return;
         }
        }        
       }
      }
   
    for(int i=0;i<nproc;i++)
   { 
    for (int s=0;s<6;s++)
      
     if (*input().getNp(s)==i) continue;
    
     ncells=selectorContainer().ProcCell(i);
     
     for (int n=0; n<ncells;n++)
      {
      MPI_Win_lock(MPI_LOCK_SHARED,i,0,*cellwin_);
        MPI_Get(&buf1,1,MPI_INT,i,n,1,MPI_INT,*cellwin_);
      MPI_Win_unlock(i,*cellwin_);
        if (buf1==cell)
        {
         buf[0]=n;
         buf[1]=i;
         return;  
        }
    
     }    
   }
  error().throw_error_one("data_storage.cpp",-1,"End of searching loop with no result");
  return;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 

void DataStorage::getRMA_V_by_coord(double* buf, int name, int* coord )
{
   int disp=*RMAvF_[name]->space();

     for (int j=0;j<3;j++) 
     {
      MPI_Win_lock(MPI_LOCK_SHARED,coord[1],0,*winV_[name][0]);
      MPI_Get(&buf[j],1,MPI_DOUBLE,coord[1],coord[0]*disp+j,1,MPI_DOUBLE,*winV_[name][0]);
      MPI_Win_unlock(coord[1],*winV_[name][0]);      
     }
}

void DataStorage::getRMA_S_by_coord(double* buf, int name, int* coord )
{
   
      MPI_Win_lock(MPI_LOCK_SHARED,coord[1],0,*winS_[name]);
      MPI_Get(buf,1,MPI_DOUBLE,coord[1],coord[0],1,MPI_DOUBLE,*winS_[name]);
      MPI_Win_unlock(coord[1],*winS_[name]);      

}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * * * */

double DataStorage::getRMA_V(int name,int component,int cell)
{
  
  int disp=*RMAvF_[name]->space();
   int ncells=selectorContainer().ProcCell(comm().me());
   
     for (int n=0; n<ncells;n++)
      {
      MPI_Win_lock(MPI_LOCK_SHARED,comm().me(),0,*cellwin_);
        MPI_Get(&buf1,1,MPI_INT,comm().me(),n,1,MPI_INT,*cellwin_);
      MPI_Win_unlock(comm().me(),*cellwin_);
        if (buf1==cell)
        {
         MPI_Win_lock(MPI_LOCK_SHARED,comm().me(),0,*winV_[name][component]);
         MPI_Get(&buf2,1,MPI_DOUBLE,comm().me(),n*disp,1,MPI_DOUBLE,*winV_[name][component]);
          MPI_Win_unlock(comm().me(),*winV_[name][component]);
         return buf2;
        }
      } 
      
      for (int i=0;i<6;i++)
      { 
        
        if (*input().getNp(i)!=-1)
        {
           ncells=selectorContainer().ProcCell(*input().getNp(i));
          
         for (int n=0; n<ncells;n++)
         {
           MPI_Win_lock(MPI_LOCK_SHARED,*input().getNp(i),0,*cellwin_);
           MPI_Get(&buf1,1,MPI_INT,*input().getNp(i),n,1,MPI_INT,*cellwin_);
           MPI_Win_unlock(*input().getNp(i),*cellwin_);
         if (buf1==cell)
         {
          MPI_Win_lock(MPI_LOCK_SHARED,*input().getNp(i),0,*winV_[name][component]);
          MPI_Get(&buf2,1,MPI_DOUBLE,*input().getNp(i),n*disp,1,MPI_DOUBLE,*winV_[name][component]);
          MPI_Win_unlock(*input().getNp(i),*winV_[name][component]);
          return buf2;
         }
        }  
        
        
        
       }
      }
   
    for(int i=0;i<nproc;i++)
   { 
    for (int s=0;s<6;s++)
      if (*input().getNp(s)==i) continue;
    ncells=selectorContainer().ProcCell(i);
     
   for (int n=0; n<ncells;n++)
      {
      MPI_Win_lock(MPI_LOCK_SHARED,i,0,*cellwin_);
        MPI_Get(&buf1,1,MPI_INT,i,n,1,MPI_INT,*cellwin_);
      MPI_Win_unlock(i,*cellwin_);
        if (buf1==cell)
        {
         MPI_Win_lock(MPI_LOCK_SHARED,i,0,*winV_[name][component]);
         MPI_Get(&buf2,1,MPI_DOUBLE,i,n*disp,1,MPI_DOUBLE,*winV_[name][component]);
         MPI_Win_unlock(i,*winV_[name][component]);
         
         
         return buf2;
        }
    
     }    
   }
  
  error().throw_error_one("data_storage.cpp",-1,"End of searching loop with no result");
  return -1;
}

double DataStorage::getRMA_S(int name,int cell)
{
   
   int ncells=selectorContainer().ProcCell(comm().me());
 
   
   for (int n=0; n<ncells;n++)
      {
      MPI_Win_lock(MPI_LOCK_SHARED,comm().me(),0,*cellwin_);
        MPI_Get(&buf1,1,MPI_INT,comm().me(),n,1,MPI_INT,*cellwin_);
      MPI_Win_unlock(comm().me(),*cellwin_);
        if (buf1==cell)
        {
         MPI_Win_lock(MPI_LOCK_SHARED,comm().me(),0,*winS_[name]);
         MPI_Get(&buf2,1,MPI_DOUBLE,comm().me(),n,1,MPI_DOUBLE,*winS_[name]);
          MPI_Win_unlock(comm().me(),*winS_[name]);
         return buf2;
        }
      } 
         
   for (int i=0;i<6;i++)
      { 
        
        if (*input().getNp(i)!=-1)
        {
        
           ncells=selectorContainer().ProcCell(*input().getNp(i));
          
         for (int n=0; n<ncells;n++)
         {
           MPI_Win_lock(MPI_LOCK_SHARED,*input().getNp(i),0,*cellwin_);
           MPI_Get(&buf1,1,MPI_INT,*input().getNp(i),n,1,MPI_INT,*cellwin_);
           MPI_Win_unlock(*input().getNp(i),*cellwin_);
         if (buf1==cell)
         {
          MPI_Win_lock(MPI_LOCK_SHARED,*input().getNp(i),0,*winS_[name]);
          MPI_Get(&buf2,1,MPI_DOUBLE,*input().getNp(i),n,1,MPI_DOUBLE,*winS_[name]);
          MPI_Win_unlock(*input().getNp(i),*winS_[name]);
          return buf2;
         }
        }  
       }
      }

   
    for(int i=0;i<nproc;i++)
   { 
   for (int s=0;s<6;s++)
      if (*input().getNp(s)==i) continue;
   ncells=selectorContainer().ProcCell(i);
   for (int n=0; n<ncells;n++)
      {
      MPI_Win_lock(MPI_LOCK_SHARED,i,0,*cellwin_);
        MPI_Get(&buf1,1,MPI_INT,i,n,1,MPI_INT,*cellwin_);
      MPI_Win_unlock(i,*cellwin_);
        if (buf1==cell)
        {
         MPI_Win_lock(MPI_LOCK_SHARED,i,0,*winS_[name]);
         MPI_Get(&buf2,1,MPI_DOUBLE,i,n,1,MPI_DOUBLE,*winS_[name]);
          MPI_Win_unlock(i,*winS_[name]);
         
         
         return buf2;
        }
    
     }    
   }
  error().throw_error_one("data_storage.cpp",-1,"End of searching loop with no result");
  return -1;
}



