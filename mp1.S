
.data					# section declaration

        # Useful offset constants for accessing members of a 
        # struct mp1_blink_struct structure
        LOCATION   = 0    # w
        ON_CHAR    = 2    # b
        OFF_CHAR   = 3    # b
        ON_LENGTH  = 4    # w
        OFF_LENGTH = 6    # w
        COUNTDOWN  = 8    # w
        STATUS     = 10   # w
        NEXT       = 12   # l

        STRUCT_SIZE = 16

# Pointer to head of list (initialized to NULL)
mp1_list_head:
        .long   0

.text					# section declaration

# Export the function symbol names

.global mp1_rtc_tasklet
.global mp1_ioctl

# void mp1_poke(void);
#
# Interface: Register-based arguments (not C-style)
#    Inputs: %cl  - The byte you wish to write
#            %eax - Offset from the start of video memory that you wish
#                   to write to
#   Outputs: Text-mode video screen is written to at location %eax with
#            the byte in %cl
# Registers: Clobbers EDX
mp1_poke:
        
	movl    vmem_base_addr(,1),%edx
	movb    %cl,(%edx,%eax,1)
	ret
# void mp1 rtc tasklet (unsigned long)
#
# Interface:
#    Inputs: no       
#   Outputs:  traverse the linked list, refresh the time property of all structs,
#              and call poke if need switch the character to be displayed for certain location
# Registers: Clobbers EAX 
mp1_rtc_tasklet:
        pushl %ebp
        movl %esp,%ebp
        movl mp1_list_head,%edx
  tk_check_end:
        cmpl $0, %edx
        je tk_return
  tk_decre_countdown:
        movw COUNTDOWN(%edx),%ax
        addw $-1, %ax   
        jz tk_see_status           # go to "see status" if "count zero"mS
        movw %ax, COUNTDOWN(%edx) # just renew countdown if not "count zero"
        movl NEXT(%edx),%edx       # %edx point to the next block 
        jmp tk_check_end           
  tk_see_status:
        movw STATUS(%edx),%ax
        cmpw $0, %ax     
        je tk_current_off
        cmpw $1, %ax     
        je tk_current_on
  tk_current_off:
        # renew COUNTDOWN and flip status
        movw $1, STATUS(%edx)
        movw ON_LENGTH(%edx),%ax
        movw %ax, COUNTDOWN(%edx)
        # call mp1_poke
        xorl %eax, %eax
        movw LOCATION(%edx), %ax  # fetch "location" to ax
        addl %eax,%eax # doulbe eax
        movb ON_CHAR(%edx),%cl    # fetch "CHAR" to cl
        pushl %edx                 # caller save: edx
        call mp1_poke
        popl %edx
        movl NEXT(%edx),%edx       # %edx point to the next block 
        jmp tk_check_end
  tk_current_on:
        # renew COUNTDOWN and flip status
        movw $0, STATUS(%edx)
        movw OFF_LENGTH(%edx),%ax
        movw %ax, COUNTDOWN(%edx)
        # call mp1_poke
        xorl %eax, %eax
        movw LOCATION(%edx), %ax  # fetch "location" to ax
        addl %eax,%eax # doulbe eax
        movb OFF_CHAR(%edx),%cl    # fetch "CHAR" to cl
        pushl %edx                 # caller save: edx
        call mp1_poke
        popl %edx
        movl NEXT(%edx),%edx       # %edx point to the next block 
        jmp tk_check_end

  tk_return:
        leave
	ret
# int mp1_ioctl (unsigned long arg, unsigned long cmd);
#
# Interface:
#    Inputs: arg (l）
#            cmd (w)         
#   Outputs:  jump to one of the inctl function
# Registers: Clobbers EAX 
mp1_ioctl:
        pushl %ebp
        movl %esp,%ebp
        
        movl 12(%ebp), %eax # fetch arg cmd to %eax
        # valid cmd test
        cmpl $0,%eax  # cmd>=0
        jl ioctl_invalid_cmd
        cmpl $3,%eax  # cmd<=3
        jg ioctl_invalid_cmd
        jmp *jump_table(,%eax,4) 
 ioctl_invalid_cmd:
        movl $-1,%eax
        leave
        ret

# int mp1 ioctl add(unsigned long arg)
#
# Interface:
#    Inputs: a user-level pointer to a mp1 blink struct structure  
#   Outputs:  copy this structure to kernal, initialize its time property, add it to linked list
# Registers: Clobbers EAX EDX ECX 
mp1_ioctl_add:
        # call mp1_malloc: allocate mem for copying the struct to kernal
        pushl $16            # arg size
        call mp1_malloc     
        addl $4,%esp         # pop args 
        # %eax = to pointer
        cmpl $0,%eax
        je add_fail_mem_allowcate 
        pushl %eax           # save "to pointer" to local:-4(%ebp) 
        movl 8(%ebp), %ecx   # fetch arg "from pointer" to %ecx
        # check if "from pointer" is NULL
        cmpl $0, %ecx
        je add_fail_mem_copy
        # call mp1_copy_from_user: copy struct from user mem to kernal mem
        pushl $16            # arg size
        pushl %ecx           # arg "from pointer"
        pushl %eax           # arg "to pointer"
        call mp1_copy_from_user
        addl $12,%esp        # pop args 
        # %eax = number of bytes fail to copy
        cmpl $0,%eax
        jne add_fail_mem_copy
        # check the location is valid
        movl -4(%ebp), %eax  # load "to pointer" to %eax
        movw LOCATION(%eax), %cx  # load "location" to %cx
        cmpw $1999, %cx     # invalid "location" >1999=80*25-1
        jg add_fail_invalid_location
        cmpw $0, %cx        # invalid "location" <0
        jl add_fail_invalid_location
        pushw %cx            # save "location" to local:-6(%ebp) 
        # set the countdown field to be equal to the on length field
        movw ON_LENGTH(%eax),%cx  # fetch "ON_LENGTH" to cx
        movw %cx, COUNTDOWN(%eax)
        # set the status field to 1
        movw $1, STATUS(%eax)
        # add this item to the linked list
        movl mp1_list_head, %ecx # fetch "head's pointer" to ecx
        movl %ecx, NEXT(%eax)
        movl %eax, mp1_list_head
        # call mp1 poke
        movl %eax, %edx     # fetch "to pointer" to %edx
        xorl %eax, %eax
        movw -6(%ebp), %ax  # load "location" to ax
        addl %eax,%eax # doulbe eax
        movb ON_CHAR(%edx),%cl
        call mp1_poke
 add_success: 
        # return 0 
        xorl %eax, %eax       
        leave
        ret

 add_fail_mem_allowcate:
        # return -1
        movl $-1,%eax
        leave
        ret
 add_fail_invalid_location:
 add_fail_mem_copy:
        # free the mem by "to pointer"
        movl -4(%ebp), %eax  # load "to pointer" to %eax
        # call mp1_free: free mem for copying the struct to kernal
        pushl %eax            # argument to pointer
        call mp1_free     
        addl $4,%esp         # pop args 
        # return -1
        movl $-1,%eax
        leave
        ret 


# mp1_blink_struct* help_locate(short location);
#
# Interface:
#    Inputs: arg1: location (w)
#            
#   Outputs:  fail:    0 (cannot find the blink_struct with this location)
#             success: pointer to the finded blink_struct
# Registers: Clobbers EAX ECX EDX 
help_locate:
        pushl %ebp
        movl %esp,%ebp
        pushl %ebx
        
        xorl %ebx,%ebx
        movw 8(%ebp), %bx        # fetch arg "location" to %bx        
        movl mp1_list_head,%edx
        movl mp1_list_head,%ecx
  hp_check_end:
        cmpl $0, %edx 
        je hp_cannot_find
  hp_compare:
        xorl %eax,%eax  
        movw LOCATION(%edx),%ax
        cmpw %bx,%ax
        je   hp_find  
        movl %edx,%ecx           # save the pointer for previous struct to %ecx
        movl NEXT(%edx),%edx     # %edx point to the next block 
        jmp hp_check_end
  hp_find:
        movl %edx,%eax
        popl %ebx
        leave
        ret     
  hp_cannot_find:
        xorl %eax,%eax
        popl %ebx
        leave
        ret

# int mp1 ioctl remove(unsigned long arg)
#
# Interface:
#    Inputs: an integer location
#   Outputs:   remove the struct with the same location from the linked list, free the mem for it
#               fail: return -1
# Registers: Clobbers EAX EDX ECX 
mp1_ioctl_remove:
        movl 8(%ebp), %ecx   # fetch arg "location" to %ecx   
        # call help_locate: to get pointer for the struct
        pushw %cx           # arg: location
        call help_locate
        addl $2, %esp      #  pop arg
        cmpl $0,%eax       # %eax store the pointer for the struct
                           # %ecx store the pointer for previous struct 
        je rm_cannot_find  
        
        # discuss if the pointer is the head
        cmpl %ecx,%eax
        je rm_head_pointer
  rm_not_head:
        # modify the linked list
        movl NEXT(%eax),%edx
        movl %edx,NEXT(%ecx)
        # free the mem for the removed struct
        # call mp1_free: free mem for copying the struct to kernal
        pushl %eax            # argument  "the pointer for the struct"
        call mp1_free     
        addl $4,%esp         # pop args 
        # return 0 
        xorl %eax, %eax       
        leave
        ret  
 rm_head_pointer:
        movl NEXT(%eax),%edx
        movl %edx,mp1_list_head
        # free the mem for the removed struct
        # call mp1_free: free mem for copying the struct to kernal
        pushl %eax            # argument  "the pointer for the struct"
        call mp1_free     
        addl $4,%esp         # pop args 
        # return 0 
        xorl %eax, %eax       
        leave
        ret 

 rm_cannot_find:
        # return -1
        movl $-1,%eax
        leave
        ret        
# int mp1 ioctl find(unsigned long arg)
#
# Interface:
#    Inputs: a pointer to a mp1 blink struct in user space
#   Outputs:   copy the blink struct with the same locaiton in kernal to the user space's structure
#               fail: return -1
# Registers: Clobbers EAX EDX ECX         
mp1_ioctl_find:
        # call mp1_malloc: allocate mem for copying the struct to kernal
        pushl $16            # arg size
        call mp1_malloc     
        addl $4,%esp         # pop args 
        # %eax = to pointer
        cmpl $0,%eax
        je find_fail_mem_allowcate 
        pushl %eax           # save "to pointer" to local:-4(%ebp) 
        movl 8(%ebp), %ecx   # fetch arg "from pointer" to %ecx
        # check if "from pointer" is NULL
        cmpl $0, %ecx
        je find_fail_mem_copy
        # call mp1_copy_from_user: copy struct from user mem to kernal mem
        pushl $16            # arg size
        pushl %ecx           # arg "from pointer"
        pushl %eax           # arg "to pointer"
        call mp1_copy_from_user
        addl $12,%esp        # pop args 
        # %eax = number of bytes fail to copy
        cmpl $0,%eax
        jne find_fail_mem_copy
        # check the location is valid
        movl -4(%ebp), %eax  # load "to pointer" to %eax
        movw LOCATION(%eax), %cx  # load "location" to %cx
        cmpw $1999, %cx     # invalid "location" >1999=80*25-1
        jg find_fail_invalid_location
        cmpw $0, %cx        # invalid "location" <0
        jl find_fail_invalid_location
        pushw %cx            # save "location" to local:-6(%ebp) 
        # call help_locate: to get pointer for the structure (in kernal)
        movw  -6(%ebp),%ax
        pushw %ax           # arg: location
        call help_locate
        addl $2, %esp      # pop arg
        cmpl $0,%eax
        je find_cannot_find  
        pushl %eax          # save "pointer(in kernal)" to local:-10(%ebp)
        # copy the entire struct in kernal to userspace
        movl -10(%ebp),%ecx  # fetch "kernal pointer" to %ecx
        movl 8(%ebp), %eax   # fetch arg "user pointer" to %eax
        pushl $16            # arg size
        pushl %ecx           # arg "kernal pointer"
        pushl %eax           # arg "user pointer"
        call mp1_copy_to_user
        addl $12,%esp        # pop args  
        # %eax = number of bytes fail to copy
        cmpl $0,%eax  
        jne find_fail_mem_copy
  find_success: 
        # free the mem by "to pointer"
        movl -4(%ebp), %eax  # load "to pointer" to %eax
        # call mp1_free: free mem for copying the struct to kernal
        pushl %eax            # argument to pointer
        call mp1_free     
        addl $4,%esp         # pop args        
        # return 0 
        xorl %eax, %eax       
        leave
        ret

 find_fail_mem_allowcate:
        # return -1
        movl $-1,%eax
        leave
        ret
 find_cannot_find:
 find_fail_invalid_location:
 find_fail_mem_copy:
        # free the mem by "to pointer"
        movl -4(%ebp), %eax  # load "to pointer" to %eax
        # call mp1_free: free mem for copying the struct to kernal
        pushl %eax            # argument to pointer
        call mp1_free     
        addl $4,%esp         # pop args 
        # return -1
        movl $-1,%eax
        leave
        ret 

        leave
	ret
# int mp1 ioctl sync(unsigned long arg)
#
# Interface:
#    Inputs: two two-byte unsigned short integers packed in a long int, upper 16bit is the locaiton of first, lower 16bits for second
#   Outputs:    copy the timing information of the first to the second
#               fail: return -1
# Registers: Clobbers EAX EDX ECX  
mp1_ioctl_sync:
        movl 8(%ebp),%eax  # fetch "first|second" argument to %eax
        pushw %ax          # save "second" to local:-2(%ebp) 
        shr $16,%eax
        pushw %ax          # save "first" to local:-4(%ebp) 

        # call help_locate: to get point for "first"
        movw  -4(%ebp),%ax
        pushw %ax           # arg: location
        call help_locate
        addl $2, %esp      # pop arg
        cmpl $0,%eax
        je sy_cannot_find_first
        pushl %eax          # save "pointer to first" to local:-8(%ebp)
        # call help_locate: to get point for "second"
        movw  -2(%ebp),%ax
        pushw %ax           # arg: location
        call help_locate
        addl $2, %esp      # pop arg
        cmpl $0,%eax
        je sy_cannot_find_second
        pushl %eax         # save "pointer to second" to local:-12(%ebp)  

        # copy time_info from "first" to "second"
        movl -8(%ebp),%ecx   # fetch "pointer to first" to %ecx
        movl -12(%ebp),%edx  # fetch "pointer to second" to %edx
        movw ON_LENGTH(%ecx),%ax
        movw %ax,ON_LENGTH(%edx)
        movw OFF_LENGTH(%ecx),%ax
        movw %ax,OFF_LENGTH(%edx)
        movw COUNTDOWN(%ecx),%ax
        movw %ax,COUNTDOWN(%edx)
        movw STATUS(%ecx),%ax
        movw %ax,STATUS(%edx)
         
        # call poke to renew display
  sy_see_status:
        cmpw $0, %ax     
        je sy_current_off
        cmpw $1, %ax     
        je sy_current_on
  sy_current_off:
        # call mp1_poke
        xorl %eax, %eax
        movw LOCATION(%edx), %ax  # fetch "location" to ax
        addl %eax,%eax # doulbe eax
        movb OFF_CHAR(%edx),%cl    # fetch "CHAR" to cl
        call mp1_poke
        jmp sy_success
  sy_current_on:
        # call mp1_poke
        xorl %eax, %eax
        movw LOCATION(%edx), %ax  # fetch "location" to ax
        addl %eax,%eax # doulbe eax
        movb ON_CHAR(%edx),%cl     # fetch "CHAR" to cl
        call mp1_poke
        jmp sy_success
 sy_success:
        # return 0
        xorl %eax, %eax
        leave
        ret         
 sy_cannot_find_first:
 sy_cannot_find_second:
        # return -1
        movl $-1,%eax
        leave
        ret       
        

/* The jump table for mp1_ioctl */
jump_table:
.long mp1_ioctl_add, mp1_ioctl_remove, mp1_ioctl_find, mp1_ioctl_sync

.end
