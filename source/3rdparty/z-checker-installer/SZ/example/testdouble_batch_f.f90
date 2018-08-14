program p
        use sz
        use rw
        implicit none
        character(len=32) :: arg
        real(kind=8), dimension(:,:,:), allocatable :: grid
        real(kind=8), dimension(:,:,:), allocatable :: grid1
        integer(kind=4) :: gridsize1,gridsize2,gridsize3
        integer(kind=4) :: g1, g2, g3, g4, g5
        real(kind=8) :: res=0
        integer :: i,j,k
        integer(kind=4) :: dim_
        integer(kind=4) :: ierr, outSize !the size of the compressed stream
        INTEGER(kind=1), DIMENSION(:), allocatable :: Bytes
        INTEGER(kind=4) :: sz_mode = 0
        gridsize1 = 10 
        gridsize2 = 10 
        gridsize3 = 10 
        
        write (6,*) 'start....'
        allocate(grid(gridsize1,gridsize2,gridsize3))
        DO i=1,gridsize3
                DO j=1,gridsize2
                        DO k=1,gridsize1
                                grid(k,j,i)=i+j+k
                        END DO
                END DO
        END DO
     
        call getarg(1, arg)
        call SZ_Init(arg,ierr)
      
        call SZ_BatchAddVar('grid1', grid, sz_mode, 0.0001D0, 0.0001D0)
        call SZ_BatchAddVar('grid2', grid, sz_mode, 0.001D0, 0.001D0)
 
        call SZ_Batch_Compress(Bytes, outSize)
        call writeData(Bytes, outSize, 'test12_f.sz') 

        !call SZ_FREE_VARSET(0)
        !deallocate(grid)

        !decompress
        call SZ_Batch_Decompress(Bytes, outSize, ierr)
        call SZ_GetVarDim('grid1', dim_, g1, g2, g3, g4, g5)
        write (6,*) 'dimension of grid1 is ',dim_
        allocate(grid1(g1,g2,g3))
        call SZ_GetVarData('grid1', grid1)
        
        open(unit=10,file='grid1_f.txt')
        DO i=1,g3
                DO j=1,g2
                        DO k=1,g1
                                write (10,*) grid1(k,j,i)
                        END DO
                END DO
        END DO

        ! Free memory
        !deallocate(grid)
        !deallocate(grid1)
        deallocate(Bytes)
        !call SZ_FREE_VARSET(1)
        !call SZ_Finalize()
        write (6,*) 'done.'
end program p
