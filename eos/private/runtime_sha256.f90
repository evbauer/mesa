module runtime_sha256
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none
   private
   public :: file_sha256
   integer(int64), parameter :: mask32=int(z'ffffffff',int64)
   integer(int64), parameter :: round_constant(64)=[ &
      int(z'428a2f98',int64),int(z'71374491',int64),int(z'b5c0fbcf',int64),int(z'e9b5dba5',int64), &
      int(z'3956c25b',int64),int(z'59f111f1',int64),int(z'923f82a4',int64),int(z'ab1c5ed5',int64), &
      int(z'd807aa98',int64),int(z'12835b01',int64),int(z'243185be',int64),int(z'550c7dc3',int64), &
      int(z'72be5d74',int64),int(z'80deb1fe',int64),int(z'9bdc06a7',int64),int(z'c19bf174',int64), &
      int(z'e49b69c1',int64),int(z'efbe4786',int64),int(z'0fc19dc6',int64),int(z'240ca1cc',int64), &
      int(z'2de92c6f',int64),int(z'4a7484aa',int64),int(z'5cb0a9dc',int64),int(z'76f988da',int64), &
      int(z'983e5152',int64),int(z'a831c66d',int64),int(z'b00327c8',int64),int(z'bf597fc7',int64), &
      int(z'c6e00bf3',int64),int(z'd5a79147',int64),int(z'06ca6351',int64),int(z'14292967',int64), &
      int(z'27b70a85',int64),int(z'2e1b2138',int64),int(z'4d2c6dfc',int64),int(z'53380d13',int64), &
      int(z'650a7354',int64),int(z'766a0abb',int64),int(z'81c2c92e',int64),int(z'92722c85',int64), &
      int(z'a2bfe8a1',int64),int(z'a81a664b',int64),int(z'c24b8b70',int64),int(z'c76c51a3',int64), &
      int(z'd192e819',int64),int(z'd6990624',int64),int(z'f40e3585',int64),int(z'106aa070',int64), &
      int(z'19a4c116',int64),int(z'1e376c08',int64),int(z'2748774c',int64),int(z'34b0bcb5',int64), &
      int(z'391c0cb3',int64),int(z'4ed8aa4a',int64),int(z'5b9cca4f',int64),int(z'682e6ff3',int64), &
      int(z'748f82ee',int64),int(z'78a5636f',int64),int(z'84c87814',int64),int(z'8cc70208',int64), &
      int(z'90befffa',int64),int(z'a4506ceb',int64),int(z'bef9a3f7',int64),int(z'c67178f2',int64)]
contains
   subroutine file_sha256(path, digest, ierr)
      character(len=*), intent(in) :: path
      character(len=64), intent(out) :: digest
      integer, intent(out) :: ierr
      character(len=65536) :: buffer
      character(len=128) :: tail
      integer(int64) :: state(8), file_size, offset, bit_length
      integer :: unit, count, position, remainder, padded, byte_index, code
      digest=''
      open(newunit=unit,file=path,status='old',access='stream',form='unformatted',action='read',iostat=ierr)
      if (ierr /= 0) return
      inquire(unit=unit,size=file_size,iostat=ierr)
      if (ierr /= 0 .or. file_size < 0 .or. file_size > huge(file_size)/8) then
         close(unit)
         ierr=1
         return
      end if
      state=[int(z'6a09e667',int64),int(z'bb67ae85',int64),int(z'3c6ef372',int64), &
         int(z'a54ff53a',int64),int(z'510e527f',int64),int(z'9b05688c',int64), &
         int(z'1f83d9ab',int64),int(z'5be0cd19',int64)]
      offset=0
      remainder=0
      tail=repeat(achar(0),128)
      do while (offset < file_size)
         count=int(min(int(len(buffer),int64),file_size-offset))
         read(unit,pos=offset+1,iostat=ierr) buffer(:count)
         if (ierr /= 0) then
            close(unit)
            return
         end if
         do position=1,count-63,64
            call compress(buffer(position:position+63),state)
         end do
         remainder=mod(count,64)
         if (remainder > 0) tail(:remainder)=buffer(count-remainder+1:count)
         offset=offset+count
      end do
      close(unit)
      tail(remainder+1:remainder+1)=achar(128)
      padded=64
      if (remainder >= 56) padded=128
      bit_length=8*file_size
      do byte_index=0,7
         tail(padded-byte_index:padded-byte_index)=achar(int(iand(shiftr(bit_length,8*byte_index),255_int64)))
      end do
      call compress(tail(:64),state)
      if (padded == 128) call compress(tail(65:128),state)
      do position=1,8
         write(digest(8*position-7:8*position),'(z8.8)') state(position)
      end do
      do position=1,64
         code=iachar(digest(position:position))
         if (code >= iachar('A') .and. code <= iachar('F')) digest(position:position)=achar(code+32)
      end do
      ierr=0
   end subroutine file_sha256

   pure function rotate(word, amount) result(rotated)
      integer(int64), intent(in) :: word
      integer, intent(in) :: amount
      integer(int64) :: rotated
      rotated=ishftc(word,-amount,32)
   end function rotate

   pure subroutine compress(block, state)
      character(len=64), intent(in) :: block
      integer(int64), intent(inout) :: state(8)
      integer(int64) :: words(64), work(8), sigma0, sigma1, choice, majority, first, second
      integer :: position, byte_index
      do position=1,16
         words(position)=0
         do byte_index=1,4
            words(position)=shiftl(words(position),8)+iachar(block(4*position-4+byte_index:4*position-4+byte_index))
         end do
      end do
      do position=17,64
         sigma0=ieor(ieor(rotate(words(position-15),7),rotate(words(position-15),18)),shiftr(words(position-15),3))
         sigma1=ieor(ieor(rotate(words(position-2),17),rotate(words(position-2),19)),shiftr(words(position-2),10))
         words(position)=iand(words(position-16)+sigma0+words(position-7)+sigma1,mask32)
      end do
      work=state
      do position=1,64
         sigma1=ieor(ieor(rotate(work(5),6),rotate(work(5),11)),rotate(work(5),25))
         choice=ieor(iand(work(5),work(6)),iand(not(work(5)),work(7)))
         first=iand(work(8)+sigma1+choice+round_constant(position)+words(position),mask32)
         sigma0=ieor(ieor(rotate(work(1),2),rotate(work(1),13)),rotate(work(1),22))
         majority=ieor(ieor(iand(work(1),work(2)),iand(work(1),work(3))),iand(work(2),work(3)))
         second=iand(sigma0+majority,mask32)
         work=[iand(first+second,mask32),work(1),work(2),work(3),iand(work(4)+first,mask32), &
            work(5),work(6),work(7)]
      end do
      state=iand(state+work,mask32)
   end subroutine compress
end module runtime_sha256
