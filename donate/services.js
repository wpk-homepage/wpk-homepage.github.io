'use strict';
/**
 * Copyright 2015 Paul Koppen,
 * released under MIT licence.
 */

angular.module('services', [])


/**
 * A provider for BTC donations.
 */
.provider('donation', [
  function() {
    
    // The default donation sends 0.1 BTC to me! MWHAHAHAHA!!!!
    var donation = {
      address: '17jEGhLXD3Tz8Rkhd8QnwGBCnmXXYnQRHU',
      amount: 0.1,
      message: 'by default.',
    };
    
    /**
     * Set the destination address.
     *
     * @param address  String of the destination bitcoin address.
     */
    this.setAddress = function(address) {
      donation.address = address;
    };
    
    /**
     * Set the amount.
     *
     * @param amount   Float specifying the amount of BTC to send.
     */
    this.setAmount = function(amount) {
      donation.amount = amount;
    };
    
    /**
     * Attach a message to the transaction.
     *
     * @param message  String to add to the transaction.
     */
    this.setMessage = function(message) {
      donation.message = message;
    };
    
    /**
     * AngularJS getter for the donation object.
     *
     * @return donation  Object with fields .address, .amount, and .message.
     */
    this.$get = [function() {
      return donation;
    }];
  }
]);