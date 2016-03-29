'use strict';
/**
 * A Bitcoin donation form.
 * Copyright 2015 Paul Koppen,
 * released under MIT license.
 */

angular.module('wpk.donate', ['wpk.qrcode', 'services'])


.config(['donationProvider', function(donationProvider){
  donationProvider.setAddress('17jEGhLXD3Tz8Rkhd8QnwGBCnmXXYnQRHU');
  donationProvider.setAmount(0.01);
  donationProvider.setMessage('');
}])

/**
 * Controller for the Donation app.
 * Requires configuration of the donationProvider.
 */
.controller('DonationController', ['$scope', 'donation',
  function($scope, donation) {
  
    // Bitcoin destination address.
    $scope.address  = donation.address;
    // Amount of BTC to send.
    $scope.amount   = donation.amount;
    // Optional message to attach (140 char max).
    $scope.message  = donation.message;
    
    // The above parameters converted to bitcoin URI, for use in the QR code.
    $scope.btcuri   = '';
    
    
    /**
     * Generate bitcoin URI from components.
     *
     * @param address  String of the destination address.
     * @param amount   Float providing the amount of BTC to send.
     * @param message  String to add a message to the transaction.
     */
    function bitcoinurl(address, amount, message) {
      return 'bitcoin:' + address + '?amount=' + amount + '&message=' + message;
    }
    
    /**
     * Update (refresh) the btcuri scope variable.
     */
    function update_btcuri() {
      $scope.btcuri = bitcoinurl($scope.address, $scope.amount, $scope.message);
    }
    
    /*
     * Watch for changes in amount and message to update the bitcoin URI.
     */
    $scope.$watch('amount', update_btcuri);
    $scope.$watch('message', update_btcuri);
  }
]);