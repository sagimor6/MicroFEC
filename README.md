A small FEC (Forward Error Correction) library, which uses minimal amount of memory.  
  
suppose:  
N = number of actual info packets you want to send in a block.  
K = number of redundancy packets you need to send in a block.  
L = fixed packet length.  
k = number of info packets lost in the transmission.  
  
structs:  
fec_inv_cache_t - serves as cache for inverses. O(N + K) space.  
fec_tx_state_t - O(N) space.  
fec_rx_state_t - O(N + K) space.  
  
functions:  
fec_inv_cache_init - O(N + K)  
fec_tx_init - O(1)  
fec_tx_add_info_pak - O(1)  
fec_tx_get_redundancy_pak - O(N\*L)  
fec_rx_init - O(1)  
fec_rx_add_pak - O(1)  
fec_rx_fill_missing_paks - O(k\*N\*L)  
fec_rx_get_info_paks - O(1)  
fec_rx_reset - O(N + K)  
  
currently we support:  
  
N + K <= 2\*\*16 + 1  
L needs to be even  
N != 0 and K != 0  
at least N packets needs to be received.  
  
TODO: make pdf to explain the math.  

